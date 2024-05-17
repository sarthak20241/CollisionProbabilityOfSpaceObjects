from OUP import *
import astropy.units as u 
import datetime 
from OUP.utils import create_poliastro_object
from time import time
import pandas as pd
import numpy as np

class MC():
    def __init__(self, verbose=False, store=False):
        self.store = store
        self.verbose = verbose

    def propagate(self, sample_points, time_range, norad_ids=None, time_unit=u.second):
        if not isinstance(sample_points, dict):
            return self._propagate_norad(sample_points, time_range, time_unit)
        if not isinstance(norad_ids, list) and norad_ids != None:
            return self._propagate_norad(sample_points[norad_ids], time_range, time_unit)
        norad_ids = sample_points.keys()

        norad_propagation = {}
        for norad in norad_ids:
            if self.verbose:
                print("MC " + norad)
            # TODO: check if norad is present in dict or not
            # TODO:  add store option
            norad_propagation[norad] = self._propagate_norad(sample_points[norad], time_range, time_unit)

        return norad_propagation                

    def _propagate_norad(self, sample_points, time_range, time_unit):
        orbit_samples = []
        for index, row in sample_points.iterrows():
            orbit_samples.append(create_poliastro_object(row[["x", "y", "z"]], row[["vx", "vy", "vz"]]))

        current = time_range[0]
        total_time = 0
        epoch_times = {}

        total_propagation = []
        mean_propagation = []
        for current in time_range:
            epoch_time = 0
            propagated_rv = []
            for i in range(len(orbit_samples)):
                temp = time()
                orbit_temp = orbit_samples[i].propagate(current*time_unit)
                epoch_time += time()-temp

                
                r = orbit_temp.r.value
                v = orbit_temp.v.value
                total_propagation.append([current, r[0], r[1], r[2], v[0], v[1], v[2]])
                propagated_rv.append([r[0], r[1], r[2], v[0], v[1], v[2]])
                

            if self.verbose:
                logger.info('MC Propagation of sample points done at %s second with time taken %s', str(current), str(epoch_time))
        
            mean_rv = np.mean(propagated_rv, axis=0)
            mean_propagation.append([current] + [mean_rv[i] for i in range(len(mean_rv))])

            total_time += epoch_time
            epoch_times[current] = epoch_time
            # current += time_delta
        
        propagated_sample_points_df = pd.DataFrame(total_propagation, columns=["timestep", "x", "y", "z", "vx", "vy", "vz"])
        propagated_mean_points_df = pd.DataFrame(mean_propagation, columns=["timestep", "x", "y", "z", "vx", "vy", "vz"])
        return propagated_sample_points_df, propagated_mean_points_df, {"total_execution_time": total_time, "execution_time_for_step":epoch_times}

        