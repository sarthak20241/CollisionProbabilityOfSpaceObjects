import sys
import os

# Add the parent directory to the sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from OUP import ESPT
import pandas as pd
import numpy as np

def seed_sample_points(mean_position, mean_velocity, std_dev_position, std_dev_velocity, num_points):
	# Generate clustered data
	np.random.seed(42)  # for reproducibility
	position_values = np.random.normal(mean_position, std_dev_position, (num_points, 3))
	velocity_values = np.random.normal(mean_velocity, std_dev_velocity, (num_points, 3))

	# Create DataFrame
	data = {
		'x': position_values[:, 0],
		'y': position_values[:, 1],
		'z': position_values[:, 2],
		'vx': velocity_values[:, 0],
		'vy': velocity_values[:, 1],
		'vz': velocity_values[:, 2]
	}

	return pd.DataFrame(data)

def run_test():
	# Number of points
	num_points = 1000

	# Mean values for clustering
	mean_position = [6328.997999, -3043.598143, -12.196120]
	mean_velocity = [1.023111, 2.098719, 7.158965]

	# Standard deviation for clustering
	std_dev_position = 10
	std_dev_velocity = 20

	sample_points = seed_sample_points(mean_position, mean_velocity, std_dev_position, std_dev_velocity, num_points)
	espt_obj = ESPT.ESPT(verbose=True)

	time_interval = [i for i in range(10)]
	df, mean_df, time = espt_obj.propagate(sample_points, time_interval)

	# print(df)
	# print(mean_df)
	# print(time)

	return df, mean_df, time


run_test()