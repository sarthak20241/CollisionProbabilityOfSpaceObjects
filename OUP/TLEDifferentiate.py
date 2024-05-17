from sgp4.api import Satrec
from sgp4.conveniences import sat_epoch_datetime
from astropy.time import Time
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, TEME, GCRS
import astropy.units as u
import numpy as np
import pandas as pd
import os
from OUP import *

class TLEDifferentiate():
    def __init__(self, verbose=False):
        self.dir = "rv_data"
        self.verbose = verbose
        use_data_path(self.dir)

    def differentiate(self, data=None, norads=None, time_range=None, timestamp=None, reference_frame="gcrs", save=False):
        if data == None and (norads == None or time_range == None):
            raise "Provide TLE data or norads and time range"
        
        if norads == None:
            norads = data.get_norads()
        elif data != None:
            norad_set = set(data.get_norads())
            for norad in norads:
                if norad not in norad_set:
                    raise "Provided norad not in data"
                

        if time_range == None:
            time_range = data.get_timerange()
        elif data != None:
            if not (time_range[0] >= data.get_timerange()[0] and time_range[1] <= data.get_timerange()[1]):
                raise "Provided time range not in data"
        
        if timestamp == None:
            timestamp = time_range[0]

        if data == None:
            means_r = {}
            means_v = {}
            covariance_r = {}
            covariance_v = {}
            rvs = {}

            for norad in norads:
                output = self._load_rv(norad, time_range, timestamp)
                if output == None:
                    raise "Data not found"
                means_r[norad] = output[0]
                means_v[norad] = output[1]
                covariance_r[norad] = output[2]
                covariance_v[norad] = output[3]
                rvs[norad] = output[4]

            return (means_r, means_v, covariance_r, covariance_v, rvs)
        
        

        return self._differentiate(data, norads, time_range, timestamp, reference_frame, save)

    def _differentiate(self, data, norads, time_range, timestamp, reference_frame, save):        
        norad_rv = {}
        norad_mean_r = {}
        norad_mean_v = {}
        norad_covariance_r = {}
        norad_covariance_v = {}
        timestamp_jd = Time(timestamp.strftime(date_format))
        for norad in norads:
            output = self._load_rv(norad, time_range, timestamp)
            if output != None:
                norad_mean_r[norad] = output[0]
                norad_mean_v[norad] = output[1]
                norad_covariance_r[norad] = output[2]
                norad_covariance_v[norad] = output[3]
                norad_rv[norad] = output[4]

            rvs = []
            for tle in data[norad]:
                parsed_tle = Satrec.twoline2rv(tle[0], tle[1])
                epoch = sat_epoch_datetime(parsed_tle).date()
                if  epoch >= time_range[0] and epoch <= time_range[1]:
                    e, r, v = parsed_tle.sgp4(timestamp_jd.jd1, timestamp_jd.jd2)

                    coords= CartesianRepresentation(r << u.km, xyz_axis=-1, differentials=CartesianDifferential(v << (u.km / u.s),xyz_axis=-1))
                    if reference_frame == "gcrs":
                        coords = TEME(coords, obstime=timestamp_jd).transform_to(GCRS(obstime=timestamp_jd)).cartesian
                    

                    rv = np.concatenate((coords.xyz.value, coords.differentials['s'].d_xyz.value))
                    rvs.append(rv.tolist())
            norad_rv[norad] = pd.DataFrame(rvs, columns=['x', 'y', 'z', 'vx', 'vy', 'vz'])
            norad_mean_r[norad] = np.mean(norad_rv[norad][["x", "y", "z"]].values, axis=0)
            norad_covariance_r[norad] = np.cov(norad_rv[norad][["x", "y", "z"]].values.T)
            norad_mean_v[norad] = np.mean(norad_rv[norad][["vx", "vy", "vz"]].values, axis=0)
            norad_covariance_v[norad] = np.cov(norad_rv[norad][["vx", "vy", "vz"]].values.T)
        
        if save:
            self._save_rvs(norad_rv, time_range, timestamp)

        return (norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_rv)

    def _load_rv(self, norad, time_range, timestamp):
        file_path = use_data_path(os.path.join(self.dir, self._rv_file_format(norad, time_range[0].strftime(date_format), time_range[1].strftime(date_format), timestamp.strftime(date_format))))
        if not os.path.exists(file_path):
            return None
    
        norad_rv = pd.read_csv(file_path)
        norad_mean_r = np.mean(norad_rv[["x", "y", "z"]].values, axis=0)
        norad_covariance_r = np.cov(norad_rv[["x", "y", "z"]].values.T)
        norad_mean_v = np.mean(norad_rv[["vx", "vy", "vz"]].values, axis=0)
        norad_covariance_v = np.cov(norad_rv[["vx", "vy", "vz"]].values.T)
        return (norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_rv)
    
    def _save_rvs(self, rv, time_range, timestamp):
        for (norad, rv) in rv.items():
            file_path = use_data_path(os.path.join(self.dir, self._rv_file_format(norad, time_range[0].strftime(date_format), time_range[1].strftime(date_format), timestamp.strftime(date_format))))
            rv.to_csv(file_path)
        


    def _rv_file_format(self, norad, start_date, end_date, timestamp):
        return "rv_%s_%s_%s_%s.csv" % (norad, start_date, end_date, timestamp)




