from ftplib import FTP
from datetime import datetime, timedelta
from pathlib import Path
from ionotec.download.extract import extract_and_remove_gz

def stations_available_ibge(start, end):

    ftp = FTP("geoftp.ibge.gov.br", timeout=60)
    ftp.login()

    list_stations = set()

    d = start

    while d <= end:

        print (d)

        year = d.year
        doy = d.timetuple().tm_yday

        remote_dir = (
            "/informacoes_sobre_posicionamento_geodesico/"
            f"rbmc/dados_RINEX3/{year}/{doy:03d}"
        )

        try:
            ftp.cwd(remote_dir)
            files = ftp.nlst()

        except Exception as e:
            print(f"ERROR {d}: {e}")
            d += timedelta(days=1)
            continue

        for fname in files:

            fname = fname.split("/")[-1]

            # Only RINEX 3 files
            if "_R_" not in fname:
                continue

            # Station ID
            station = fname[:9].upper()

            list_stations.add(station)

        d += timedelta(days=1)

    ftp.quit()

    return sorted(list_stations)


def download_ibge(station, start, end, output_dir="."):

    station = station.upper()
    site = station[:4]
    print (output_dir)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    ftp = FTP("geoftp.ibge.gov.br", timeout=60)
    ftp.login()

    d = start

    while d <= end:

        year = d.year
        doy = d.timetuple().tm_yday

        remote_dir = (
            f"/informacoes_sobre_posicionamento_geodesico/"
            f"rbmc/dados_RINEX3/{year}/{doy:03d}"
        )

        try:
            ftp.cwd(remote_dir)
            files = ftp.nlst()

        except Exception as e:
            print(f"ERROR accessing {remote_dir}: {e}")
            d += timedelta(days=1)
            continue

        found = False

        for fname in files:

            fname = fname.split("/")[-1]

            # Match the complete station identifier
            if not fname.upper().startswith(station):
                continue

            # We want compressed RINEX files
            if not fname.endswith(".gz"):
                continue

            found = True

            gz_path = output_dir / fname
            crx_path = gz_path.with_suffix("")

            # Already extracted
            if crx_path.exists():
                print(f"SKIP  {crx_path.name}")
                list_files.append(crx_path)
                continue

            #print(f"GET   {fname}")

            try:
                with open(gz_path, "wb") as f:
                    ftp.retrbinary(
                        f"RETR {fname}",
                        f.write
                    )

                # Extract .gz -> .crx
                crx_path = extract_and_remove_gz(gz_path)

                list_files.append(crx_path)

                #print(f"      -> {crx_path.name}")

            except Exception as e:

                if gz_path.exists():
                    gz_path.unlink()

                print(f"ERROR {fname}: {e}")

        if not found:
            print(f"NONE  {station} on {d} (DOY {doy:03d})")

        d += timedelta(days=1)

    ftp.quit()

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files
