from OUP import *
import numpy as np
import pandas as  pd
from datetime import datetime, timedelta
from scipy.stats import multivariate_normal
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import astropy.units as u 
from astropy.time import Time


date_format = "%Y-%m-%d"
class ParticleSampling:
    # we will be given arrays as follows
    # has to generate sample_count amount of samples
    # Update it as po

    def __init__(self, save=False):
        self.save = save
        self.dir = 'particle_sampling_data'

    # TODO: for multiple distribution types
    def generate_samples(
        self,
        norad_mean_r,
        norad_mean_v,
        norad_covariance_r,
        norad_covariance_v,
        norad_ids,
        sample_count,
        distribution_type="multivariate_normal",
        timestamp=''
    ):
        norad_sample_points = {}

        for id in norad_ids:
            mean_r = norad_mean_r[id]
            mean_v = norad_mean_v[id]
            covariance_r = norad_covariance_r[id]
            covariance_v = norad_covariance_v[id]
            
            # options for multiple types of distribution
            # TODO: add other distribution options here
            pdf_r = None
            pdf_v = None
            if distribution_type=="multivariate_normal":
                pdf_r = multivariate_normal(mean_r, covariance_r)
                pdf_v = multivariate_normal(mean_v, covariance_v)
            

            # Takes samples from multivariate normal distribution
            samples_r = pdf_r.rvs(size=sample_count)
            samples_v = pdf_v.rvs(size=sample_count)
            
            sample_point = []
            
            for i in range(sample_count):
                point = []
                r = samples_r[i, 0:3]*u.km
                point.extend(r.value)
                v = samples_v[i, 0:3]*u.km/u.s
                point.extend(v.value)
                sample_point.append(point)
            
            sample_point_df = pd.DataFrame(sample_point, columns=['x','y','z','vx','vy','vz'])
            norad_sample_points[id] = sample_point_df

        if self.save:
            self._save_norad_sample_points(norad_sample_points, timestamp)
        
        return norad_sample_points
    
    def _save_norad_sample_points(self, norad_sample_points, timestamp):
        for (id, sample_point_df) in norad_sample_points.items():
            timestamp = '_' + timestamp if timestamp else ''
            file_name = f'sample_points_{id}{timestamp}.csv'
            path = use_data_path(self.dir + '/' + file_name)
            sample_point_df.to_csv(path)
