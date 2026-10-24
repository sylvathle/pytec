from pathlib import Path
from datetime import timedelta

import time
import subprocess
import requests

# RAMSAC (IGN Argentina) daily RINEX 2 Hatanaka files, Unix-compressed (.Z):
#   https://ramsac.ign.gob.ar/api/v1/rinex/download/5/fosa0010.24d.Z
# The server does not send its intermediate certificate (Let's Encrypt YR1),
# so requests needs a CA bundle that includes it. Build it once, e.g.:
#   curl -sS -o yr1.der http://yr1.i.lencr.org/
#   openssl x509 -inform DER -in yr1.der -out yr1.pem
#   cat /etc/ssl/certs/ca-certificates.crt yr1.pem > ~/.cache/ionotec/ramsac-bundle.pem
BASE_URL = "https://ramsac.ign.gob.ar/api/v1/rinex/download"
CA_BUNDLE = Path.home() / ".cache" / "ionotec" / "ramsac-bundle.pem"

# Selector that appears in the download URL (5 works for the daily files;
# 1 gives a 500 error and unknown values give 403). Tried in this order.
IDS = (5,)

# Status codes that mean "this file does not exist"
_NONE_CODES = (403, 404, 500)

# Station codes of the RAMSAC status page (lower case)
STATIONS = [s.lower() for s in """
25MA 3ARO ABRA AGGO ALTA ALUM AUTF AZUL BATE BCAR BCH1 BOLS BORC BSEN CAEP CAQI
CATA CFAG CGSJ CHAC CHIM CHLT CHMA CHOY CORD COYQ CSJ1 CSLO CSOM DINO DORE EBY1
ECAR EPGZ EPSF ESQU EZZA FEDE FMAT FOSA FUTF GGUA GRCA GRLS GUAY GVIL IGM1 ITAI
JBAL JCBC JCHL JVCH JVGO LARJ LHCL LMHS LNQM LPGS LROS MBIO MECO MGLO MGNO MGUE
MPL2 MZAC MZAE MZAL MZAR MZAU MZGA MZRF NEQN NESA NOYA NUGA OAFA OBRA OJOA ORAN
PACO PATA PBCA PDE3 PEBA PEJ1 PELU PIAG PIBV PORT PRCV PRNA PUMA PWIL PWRO PY01
RAF1 RCRO RDCM RECO RIO2 RIO4 RODE RSAL RSAU RSCL RUFI RWSN SAIS SAJA SANT SARM
SBAL SICO SJSF SL01 SMAN SMDM SPRZ SRLP SUAR SUR1 SVIC TATU TAVA TERO TGTA TILC
TOLH TOSF TRN1 TUC1 TUCU UCOR UNPA UNRO UNSA UNSJ UYBU UYCO UYPA UYSA UYSO UYTD
VBCA VCON VCOP VDMA VIMA VIME VIRA VLAN YCBA YEMA ZPLA
""".split()]


def _new_session():
    s = requests.Session()
    s.verify = str(CA_BUNDLE) if CA_BUNDLE is not None and Path(CA_BUNDLE).exists() else True
    return s


def _filename(station, d):
    """ssssDDD0.YYd.Z  (daily file, day of year, 2-digit year)"""
    return f"{station}{d.timetuple().tm_yday:03d}0.{d.year % 100:02d}d.Z"


def _call(session, method, url, max_retries=4, **kwargs):
    """Request with retries on connection problems. Returns the Response."""
    for attempt in range(max_retries):
        try:
            return session.request(method, url, timeout=(10, 120), **kwargs)

        except requests.exceptions.SSLError as e:
            raise RuntimeError(
                "TLS certificate error: ramsac.ign.gob.ar does not send its "
                "intermediate certificate. Build the CA bundle described at "
                f"the top of argentina.py (expected at {CA_BUNDLE}). {e}"
            )

        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.Timeout,
        ) as e:
            if attempt == max_retries - 1:
                raise

            wait = 2 * (attempt + 1)
            print(f"WARNING {url}: {e} -> retry in {wait} s")
            time.sleep(wait)


def _exists(session, station, d):
    """The selector ID of the first URL that has the day's file, or None.
    Uses HEAD (no download); falls back to a streamed GET if HEAD is refused."""
    fname = _filename(station, d)

    for sid in IDS:
        url = f"{BASE_URL}/{sid}/{fname}"

        r = _call(session, "HEAD", url, allow_redirects=True)

        if r.status_code in (405, 501):
            r = _call(session, "GET", url, stream=True)
            r.close()

        if r.status_code == 200:
            return sid

        if r.status_code not in _NONE_CODES:
            r.raise_for_status()

    return None


def _uncompress(z_path):
    """Extract a Unix-compressed .Z file (gzip -dc also reads .Z); remove
    the .Z afterwards. Returns the extracted Path."""
    z_path = Path(z_path)
    out_path = z_path.with_suffix("")

    for cmd in (["gzip", "-dc", str(z_path)], ["uncompress", "-c", str(z_path)]):
        try:
            with open(out_path, "wb") as f:
                subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)

            z_path.unlink()
            return out_path

        except FileNotFoundError:
            continue

        except subprocess.CalledProcessError as e:
            out_path.unlink(missing_ok=True)
            raise RuntimeError(f"{cmd[0]} failed on {z_path.name}: {e.stderr!r}")

    out_path.unlink(missing_ok=True)
    raise RuntimeError("neither gzip nor uncompress was found in the PATH")


def stations_available_argentina(start, end, stations=None):
    """
    Return the union of RAMSAC stations having a daily RINEX file during
    the requested period (4-character lower-case names, like the other
    sources).

    There is no listing endpoint that I know of, so every candidate station
    is checked for every day with a HEAD request (no download). With the
    ~140 stations of the network that is ~140 requests per day: pass
    `stations=[...]` to restrict the candidates.

    Archive:
        https://ramsac.ign.gob.ar/api/v1/rinex/download/5/ssssDDD0.YYd.Z
    """

    candidates = [s.lower() for s in (stations or STATIONS)]

    list_stations = set()

    with _new_session() as session:

        d = start

        while d <= end:

            for st in candidates:

                if st in list_stations:
                    continue

                try:
                    if _exists(session, st, d) is not None:
                        list_stations.add(st)

                except RuntimeError:
                    raise

                except Exception as e:
                    print(f"ERROR {st} {d}: {e}")

            d += timedelta(days=1)

    return sorted(list_stations)


def download_argentina(station, start, end, output_dir="."):
    """
    Download RAMSAC (IGN Argentina) daily RINEX 2 Hatanaka files.

    Archive structure:
        https://ramsac.ign.gob.ar/api/v1/rinex/download/5/ssssDDD0.YYd.Z

    station: 4-character name (case-insensitive; 9-character ids are cut
    to their first 4 characters).

    The downloaded .Z files are extracted with gzip (or uncompress) and
    removed; the files are kept as Hatanaka RINEX 2 (.YYd).

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """

    station = station.lower()

    if len(station) == 9:
        station = station[:4]

    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    with _new_session() as session:

        d = start

        while d <= end:

            doy = d.timetuple().tm_yday
            fname = _filename(station, d)

            z_path = output_dir / fname
            rinex_path = z_path.with_suffix("")

            # Already extracted
            if rinex_path.exists():
                print(f"SKIP  {rinex_path.name}")
                list_files.append(rinex_path)
                d += timedelta(days=1)
                continue

            try:
                done = False

                for sid in IDS:

                    print(f"GET   {fname}")

                    r = _call(
                        session, "GET", f"{BASE_URL}/{sid}/{fname}", stream=True
                    )

                    with r:
                        if r.status_code in _NONE_CODES:
                            continue

                        r.raise_for_status()

                        with open(z_path, "wb") as f:
                            for chunk in r.iter_content(chunk_size=1024 * 1024):
                                if chunk:
                                    f.write(chunk)

                    list_files.append(_uncompress(z_path))
                    done = True
                    break

                if not done:
                    print(f"NONE  {station} on {d} (DOY {doy:03d})")

            except RuntimeError:
                raise

            except Exception as e:

                z_path.unlink(missing_ok=True)
                rinex_path.unlink(missing_ok=True)

                print(f"ERROR {fname}: {e}")

            d += timedelta(days=1)

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files
