from spacetrack import SpaceTrackClient
import spacetrack.operators as spacetrackop
import os
import json
import datetime
from dateutil.parser import parse
from OUP.utils import date_range_list
from OUP import *

class TLEData():
    def __init__(self, data, norads, time_range):
        self.data = data
        self.time_range = time_range
        self.norads = norads

    def get_timerange(self):
        return self.time_range
    
    def get_norads(self):
        return self.norads
    
    def __getitem__(self, norad):
        return self.data[norad]
    
    def __str__(self):
        return self.data.__str__()
    
    def __repr__(self):
        return self.data.__repr__()
    
    

class TLEFetch():
    def __init__(self, username, password, verbose=False):
        self.username = username
        self.password = password
        self.dir = "tle_data"
        use_data_path(self.dir)
        self.satcat_file = "satcat.json"
        self.satcat_date_file = "satcat_date"
        self.verbose = verbose

        self.spacetrack_client = SpaceTrackClient(self.username, self.password)
    
    def get_data(self, norads=None, time_range=None):
        if norads == None:
            norads = self.get_satcat()
        if time_range == None:
            time_range = (datetime.date.today()-datetime.timedelta(weeks=1), datetime.date.today())
        
        data, complete = self._get_data_file(norads, time_range)
        if complete:
            return data
        
        for norad in norads:
            for date in date_range_list(time_range):
                file_path = use_data_path(os.path.join(self.dir, self._tle_file_format(norad, date)))
                file = open(file_path, "w")
                file.close()

        data = self._get_data_spacetrack(norads, time_range).splitlines()
        return_data = {}
        for i in range(0, len(data), 2):
            norad = self._extract_norad_from_tle(data[i])
            date = self._extract_date_from_tle(data[i])
            file_path = use_data_path(os.path.join(self.dir, self._tle_file_format(norad, date)))
            file = open(file_path, "a")
            file.write(data[i] + '\n')
            file.write(data[i+1] + '\n')
            file.close()

            if not norad in return_data:
                return_data[norad] = []
            return_data[norad].append((data[i], data[i+1]))
        return TLEData(return_data, list(return_data.keys()), time_range)

    def _get_data_file(self, norads, time_range):
        return_data = {}
        complete = True
        
        for norad in norads:
            li = []
            for date in date_range_list(time_range):
                data_file_path = use_data_path(os.path.join(self.dir, self._tle_file_format(norad, date)))
                if not os.path.exists(data_file_path):
                    complete = False
                    continue
                file = open(data_file_path, "r")

                data = file.read().splitlines()
                for i in range(0, len(data), 2):
                    li.append((data[i], data[i+1]))
                file.close()
            if len(li) > 0:
                return_data[norad] = li

        return (TLEData(return_data, list(return_data.keys()), time_range), complete)
    
    def _get_data_spacetrack(self, norads, time_range):
        (start, end) = time_range
        start = start.strftime(date_format)
        end = end.strftime(date_format)

        return self.spacetrack_client.tle(norad_cat_id=norads, epoch=spacetrackop.inclusive_range(start, end), format="tle")
                    
    
    def get_satcat(self):
        data_file_path = use_data_path(os.path.join(self.dir, self.satcat_file))
        date_file_path = use_data_path(os.path.join(self.dir, self.satcat_date_file))
        if os.path.exists(data_file_path):
            date_file = open(date_file_path)
            date = parse(date_file.readline())
            date_file.close()
            if datetime.datetime.now()-date <= datetime.timedelta(weeks=1):
                data_file = open(data_file_path, 'r')
                data_str = data_file.read()
                data = json.loads(data_str)
                data = [sat["NORAD_CAT_ID"] for sat in data]
                data_file.close()
                return data

            

        data = self.spacetrack_client.satcat()
        data_file = open(data_file_path, "w")
        data = [item for item in data if item["CURRENT"] == "Y"]
        data_file.write(json.dumps(data))
        data_file.close()

        date_file = open(date_file_path, "w")
        date_file.write(datetime.datetime.now().isoformat())
        date_file.close()
        data = [sat["NORAD_CAT_ID"] for sat in data]
        return data

    def _tle_file_format(self, norad_id, date):
        return ("tle_%s_%s.txt" % (norad_id, date))
    
    # def _date_range_list(self, time_range):
    #     (start, end) = time_range
    #     return [(start+datetime.timedelta(days=day)).strftime(self.date_format) for day in range((end-start).days+1)]

    def _extract_date_from_tle(self, line):
        epoch = line.split()[3]
        year = int(epoch[:2])
        if year >= 57:
            year += 1900
        else:
            year += 2000
        days = int(float(epoch[2:]))
        return (datetime.date(year, 1, 1) + datetime.timedelta(days=days-1)).strftime(date_format)

    def _extract_norad_from_tle(self, line):
        return line.split()[1][:-1]
