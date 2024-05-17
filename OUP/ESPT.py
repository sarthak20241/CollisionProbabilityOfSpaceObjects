from OUP import *
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from scipy.linalg import expm
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import astropy.units as u 
import time
from OUP.utils import create_poliastro_object


class ESPT:
	# TODO: Allow users to input their own dynamics functions
	# TODO: Ask whether the delta in d_delta_f is the full delta (containing both position and velocity)

	def __init__(self, verbose=False, save=False):
		self.verbose = verbose
		self.dir = 'ESPT_data'
		self.save = save

	def _gravitational_matrix(self, r):
		"""
		Defines state dynamics of a satellite around earth
		"""

		mu = 398600.4418  # Gravitational parameter for Earth in km^3/s^2

		# Identity and null matrices
		I = np.identity(3)
		zero_matrix = np.zeros((3, 3))

		# Compute the gravitational term
		r_magnitude = np.linalg.norm(r)
		gravitational_term = -(mu / (r_magnitude**3)) * I

		# Construct the matrix
		matrix_top = np.hstack([zero_matrix, I])
		matrix_bottom = np.hstack([gravitational_term, zero_matrix])
		result_matrix = np.vstack([matrix_top, matrix_bottom])

		return result_matrix

	def _jacob_matrix(self, jac_function, point):
		return jac_function(point)

	def _exp_matrix(self, _jacob_matrix, time):
		return expm(_jacob_matrix*time)

	def _d_delta_F(self, _exp_matrix, delta):
		return np.dot(_exp_matrix, delta)

	def _n1(self, F_at_mean, mean, delta, time):
		J = self._jacob_matrix(self._gravitational_matrix, mean)
		D = self._d_delta_F(self._exp_matrix(J, time), delta)
		
		return F_at_mean + D

	def _n2(self, F_at_mean, mean, delta, delta_r, time):
		J1 = self._jacob_matrix(self._gravitational_matrix, mean)
		J2 = self._jacob_matrix(self._gravitational_matrix, mean + delta_r)
		D1 = self._d_delta_F(self._exp_matrix(J1, time), delta)
		D2 = self._d_delta_F(self._exp_matrix(J2, time), delta)

		return F_at_mean + D1 + D2

	def _n3(self, F_at_mean, mean, delta, delta_r, time):
		J1 = self._jacob_matrix(self._gravitational_matrix, mean)
		J2 = self._jacob_matrix(self._gravitational_matrix, mean + delta_r)
		J3 = self._jacob_matrix(self._gravitational_matrix, mean + (2*delta_r))
		D1 = self._d_delta_F(self._exp_matrix(J1, time), delta)
		D2 = self._d_delta_F(self._exp_matrix(J2, time), delta)
		D3 = self._d_delta_F(self._exp_matrix(J3, time), delta)
		return F_at_mean + D1 + D2 + D3

	def espt_ith_sample_state(self, mean, F_at_mean, time, delta):
		half_delta = delta / 2
		half_delta_r = half_delta[0:3]

		term1 = 2*self._n2(F_at_mean, mean, half_delta, half_delta_r, time)
		term2 = self._n1(F_at_mean, mean, delta, time)

		return term1 - term2

	def aespt_ith_sample_state(self, mean, F_at_mean, time, delta):
		half_delta = delta / 2
		half_delta_r = half_delta[0:3]
		third_delta = delta / 3
		third_delta_r = third_delta[0:3]
		
		term1 = 0.5*self._n1(F_at_mean, mean, delta, time)
		term2 = -4*self._n2(F_at_mean, mean, half_delta, half_delta_r, time)
		term3 = 4.5*self._n3(F_at_mean, mean, third_delta, third_delta_r, time)

		return term1 + term2 + term3

	def _orbit_propagator(self , poliastro_object, time):
		return poliastro_object.propagate(time*u.second) 

	def _apply_operation(self, row, mean, propagated_mean_point, time, propagation_type):
		if propagation_type=='espt':
			return self.espt_ith_sample_state(mean, propagated_mean_point, time, row)
		else:
			return self.aespt_ith_sample_state(mean, propagated_mean_point, time, row)
		 

	def propagate(self, sample_points, time_interval, propagation_type='espt', norad_id=''):
		# TODO: fix time interval, enable user to set the step
		# TODO: also include what if user just wants to give a date_time only
		"""
		Parameters
		---------------
		sample_points: pd.dataframe
			coordinates of sample points with column headers x,y,z,vx,vy,vz
		time_interval: array
			interval for which the propagation needs to be done
		propagation_type: string
			define the type of propagation to use. ('espt'/'aespt')

		Returns
		-----------
		propagated_sample_points_df: Pandas Dataframe
			Dataframe containing position and velocity coordinates of propagated sample points
		propagated_mean_df: np.array: np.array
			propagated_mean of the sample points
		execution_time: float
			Time it takes for the propagation to be executed
		"""

		# Find r and v
		r = sample_points.loc[:, ['x','y','z']]
		v = sample_points.loc[:, ['vx','vy','vz']]
		
		# Calculate mean values
		mean_r = r.mean()
		mean_v = v.mean()
		
		# Create poliastro object
		mean_orbit_object = create_poliastro_object(mean_r, mean_v)

		# Calculate delta
		delta_r = r.sub(mean_r)
		delta_v = v.sub(mean_v)
		delta = pd.concat([delta_r, delta_v], axis=1)
		
		propagated_sample_points_arr = []
		propagated_mean_point_arr = []

		# stores execution time in seconds
		execution_time = {
			'total_execution_time': 0,
			'execution_time_for_step': {}
		}
		
		for t in time_interval:
			start = time.time()
			propagated_poliastro_object = self._orbit_propagator(mean_orbit_object, t)
			time_taken = time.time() - start
			
			rv_vector = np.array(propagated_poliastro_object.rv())
			propagated_mean_point = np.concatenate((rv_vector[0],rv_vector[1])) 
			new_propagated_mean_point = np.insert(propagated_mean_point, 0, t)
			propagated_mean_point_arr.append(new_propagated_mean_point)

			# temp = time.time()
			propagated_sample_points = delta.apply(self._apply_operation, args=(mean_r, propagated_mean_point, t, propagation_type), axis=1)
			# time_taken += time.time()-temp

			new_propagated_sample_point = []
			for point in propagated_sample_points:
				new_point = np.insert(point, 0, t)
				new_propagated_sample_point.append(new_point)		
			
			propagated_sample_points_arr.extend(new_propagated_sample_point)
			
			execution_time['execution_time_for_step'][t] = time_taken
			execution_time['total_execution_time'] += time_taken
			
			if (self.verbose):
				logger.info('Propagation of sample points done at %s second with time taken %s', str(t), str(time_taken))
			
		propagated_sample_points_df = pd.DataFrame(propagated_sample_points_arr, columns=['timestep', 'x','y','z','vx','vy','vz'])
		propagated_mean_point_df = pd.DataFrame(propagated_mean_point_arr, columns=['timestep', 'x','y','z','vx','vy','vz'])
		
		if self.save:
			self._save_propagated_sample_points_data(propagated_sample_points_df, propagated_mean_point_df, time_interval, norad_id)
		
		return propagated_sample_points_df, propagated_mean_point_df, execution_time
	

	def _save_propagated_sample_points_data(self, sample_points, mean_points, time_interval, id):
		'''
		Save sample point dataframe and mean point dataframe in .csv file
		'''
		sample_points_file_name = f'espt_propagated_{id}_{time_interval[0]}_{time_interval[1]}.csv'
		mean_points_file_name = f'espt_mean_propagated_{id}_{time_interval[0]}_{time_interval[1]}.csv'
		
		sample_point_path = use_data_path(self.dir + '/' + sample_points_file_name)
		mean_points_path = use_data_path(self.dir + '/' + mean_points_file_name)
		
		sample_points.to_csv(sample_point_path)
		mean_points.to_csv(mean_points_path)