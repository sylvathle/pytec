
from ftplib import FTP
from pathlib import Path

import subprocess
import requests
import re

from datetime import timedelta
from ionotec.download.extract import extract_and_remove_gz


def stations_available_sugar(start, end):

    ftp = FTP("ftp.earthobservatory.sg", timeout=60)
    ftp.login()

    list_stations = set()

    for year in range(start.year, end.year + 1):

        remote_dir = f"/SugarData/{year}"

        try:
            ftp.cwd(remote_dir)
            stations = ftp.nlst()

        except Exception as e:
            print(f"ERROR {year}: {e}")
            continue

        for station in stations:

            station = station.split("/")[-1]

            # SuGAr stations are 4-character RINEX v2 codes
            if len(station) != 4:
                continue

            list_stations.add(station)

    ftp.quit()

    return sorted(list_stations)


def download_sugar(station, start, end, output_dir="."):

    station = station.lower()

    if len(station) != 4:
        raise ValueError(
            f"SuGAr station must be a 4-character RINEX v2 code: {station}"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    list_files = []

    ftp = FTP("ftp.earthobservatory.sg", timeout=60)
    ftp.login()

    d = start

    while d <= end:

        year = d.year
        doy = d.timetuple().tm_yday

        remote_dir = (
            f"/SugarData/{year}/{station}/RINEX2-HATANAKA"
        )

        fname = f"{station}{doy:03d}0.{year % 100:02d}d.gz"

        try:
            ftp.cwd(remote_dir)

        except Exception as e:
            print(f"ERROR accessing {remote_dir}: {e}")
            d += timedelta(days=1)
            continue

        gz_path = output_dir / fname
        crx_path = gz_path.with_suffix("")

        # Already extracted
        if crx_path.exists():
            print(f"SKIP  {crx_path.name}")
            list_files.append(crx_path)
            d += timedelta(days=1)
            continue

        # Download
        try:
            print(f"GET   {fname}")

            with open(gz_path, "wb") as f:
                ftp.retrbinary(
                    f"RETR {fname}",
                    f.write
                )

        except Exception as e:

            if gz_path.exists():
                gz_path.unlink()

            print(f"NONE  {station} on {d} "
                  f"(DOY {doy:03d}): {e}")

            d += timedelta(days=1)
            continue

        # Extract .gz -> .d
        try:
            crx_path = extract_and_remove_gz(gz_path)
            list_files.append(crx_path)

        except Exception as e:

            if gz_path.exists():
                gz_path.unlink()

            print(f"ERROR extracting {fname}: {e}")

        d += timedelta(days=1)

    ftp.quit()

    #print(f"Files available locally: {len(list_files)}")

    return list_files
