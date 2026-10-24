
from ftplib import FTP
from pathlib import Path

import subprocess
import requests
import re

from datetime import timedelta
from ionotec.download.extract import extract_and_remove_gz

def stations_available_cors(start, end):
    """
    Return the union of 4-character NOAA CORS stations
    having daily Hatanaka RINEX files during the requested period.

    Archive:
        https://noaa-cors-pds.s3.amazonaws.com/rinex/YYYY/DDD/

    File:
        ssssDDD0.YYd.gz
    """

    list_stations = set()

    d = start

    with requests.Session() as session:

        while d <= end:

            year = d.year
            doy = d.timetuple().tm_yday

            prefix = f"rinex/{year}/{doy:03d}/"

            continuation_token = None

            while True:

                params = {
                    "list-type": "2",
                    "prefix": prefix,
                }

                if continuation_token is not None:
                    params["continuation-token"] = continuation_token

                try:
                    r = session.get(
                        "https://noaa-cors-pds.s3.amazonaws.com/",
                        params=params,
                        timeout=60
                    )
                    r.raise_for_status()

                except Exception as e:
                    print(f"ERROR {d}: {e}")
                    break

                # Extract object keys
                keys = re.findall(
                    r"<Key>([^<]+)</Key>",
                    r.text
                )

                # Expected:
                # rinex/YYYY/DDD/ssss/ssssDDD0.YYd.gz
                pattern = re.compile(
                    rf"^rinex/{year}/{doy:03d}/"
                    rf"([A-Za-z0-9]{{4}})/"
                    rf"[A-Za-z0-9]{{4}}{doy:03d}0\.\d{{2}}d\.gz$",
                    re.IGNORECASE
                )

                for key in keys:

                    match = pattern.match(key)

                    if match:
                        list_stations.add(
                            match.group(1).lower()
                        )

                # Is there another page?
                if "<IsTruncated>true</IsTruncated>" not in r.text:
                    break

                match = re.search(
                    r"<NextContinuationToken>([^<]+)</NextContinuationToken>",
                    r.text
                )

                if not match:
                    print(
                        f"WARNING: no continuation token for {d}"
                    )
                    break

                continuation_token = match.group(1)

            d += timedelta(days=1)

    return sorted(list_stations)


def download_cors(station, start, end, output_dir="."):
    """
    Download NOAA CORS RINEX 2 daily Hatanaka files.

    Archive structure:
        https://noaa-cors-pds.s3.amazonaws.com/rinex/YYYY/DDD/ssss/

    Filename:
        ssssDDD0.YYd.gz

    Example:
        https://noaa-cors-pds.s3.amazonaws.com/rinex/2024/096/akrd/akrd0960.24d.gz

    The downloaded .gz files are extracted using
    extract_and_remove_gz().

    Returns
    -------
    list of Path
        Locally available RINEX files.
    """

    station = station.lower()

    if len(station) != 4:
        raise ValueError("station must be a 4-character station name")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    d = start

    with requests.Session() as session:

        while d <= end:

            year = d.year
            doy = d.timetuple().tm_yday
            yy = year % 100

            fname = f"{station}{doy:03d}0.{yy:02d}d.gz"

            url = (
                f"https://noaa-cors-pds.s3.amazonaws.com/"
                f"rinex/{year}/{doy:03d}/{station}/{fname}"
            )

            gz_path = output_dir / fname
            rinex_path = gz_path.with_suffix("")

            # Already extracted
            if rinex_path.exists():
                print(f"SKIP  {rinex_path.name}")
                list_files.append(rinex_path)
                d += timedelta(days=1)
                continue

            try:
                print(f"GET   {fname}")

                with session.get(
                    url,
                    stream=True,
                    timeout=120
                ) as response:

                    if response.status_code == 404:
                        print(
                            f"NONE  {station} on {d} "
                            f"(DOY {doy:03d})"
                        )
                        d += timedelta(days=1)
                        continue

                    response.raise_for_status()

                    with open(gz_path, "wb") as f:
                        for chunk in response.iter_content(
                            chunk_size=1024 * 1024
                        ):
                            if chunk:
                                f.write(chunk)

                rinex_path = extract_and_remove_gz(gz_path)
                list_files.append(rinex_path)

            except Exception as e:

                if gz_path.exists():
                    gz_path.unlink()

                print(f"ERROR {fname}: {e}")

            d += timedelta(days=1)

    print()
    print(f"Files available locally: {len(list_files)}")

    return list_files
