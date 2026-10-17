from ionotec.download import igs
from ionotec.download import ibge

from pathlib import Path

class downrinex:

    def __init__(self, start, end, output_folder='.'):

        self.date_start = start
        self.date_end = end
        self.dict_station_source = {}
        self.list_stations = []

        self.folder_deposit = output_folder

        Path(output_folder).mkdir(parents=True, exist_ok=True)


    def list_available_stations(self,source='all'):
        self.dict_station_source = {}
        if source in ['all','ibge']:
            available_stations = ibge.stations_available_ibge(self.date_start, self.date_end)
            for station in available_stations:
                self.dict_station_source[station]='ibge'
                self.list_stations.append(station)

        return self.list_stations

    def download_rinex_station(self, station):
        source = self.dict_station_source[station]
        if source  == 'ibge':
            return ibge.download_ibge(station, self.date_start, self.date_end, output_dir=self.folder_deposit)


