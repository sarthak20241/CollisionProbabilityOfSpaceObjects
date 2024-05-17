import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)


from OUP import ParticleSampling
from OUP.TLEData import TLEFetch
from OUP.TLEDifferentiate import TLEDifferentiate
from OUP.ParticleSampling import ParticleSampling
from OUP.ESPT import ESPT
from OUP.MC import MC
from OUP.Visualize import Visualize
import datetime

time_range_tle = (datetime.date(2023, 4, 16), datetime.date(2023, 4, 30))
norad_ids = ["42962"]

#TLE Fetching
fetcher_obj = TLEFetch("sanat@iiitd.ac.in", "DhgiQQvdVlOfHT49")
data = fetcher_obj.get_data(norad_ids, time_range=time_range_tle)

#TLE differentiating
diff_obj = TLEDifferentiate()
norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_rv = diff_obj.differentiate(data, save=True)

particle_obj = ParticleSampling(True)
norad_samples = particle_obj.generate_samples(norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_ids, 100)


mc_obj = MC(False, False)
espt_obj = ESPT(verbose=False, save=False)

for (norad_id, sample_points) in norad_samples.items():
	df_mc, mean_df_mc, execution_time_mc = mc_obj.propagate(sample_points, range(0,6000,600), norad_ids=norad_id)
	df_espt, mean_df_espt, execution_time_espt = espt_obj.propagate(sample_points, range(0,6000,600), norad_id=norad_id)
	df_aespt, mean_df_aespt, execution_time_aespt = espt_obj.propagate(sample_points, range(0,6000,600), norad_id=norad_id, propagation_type="aespt")
	
	Visualize.execution_time([execution_time_mc, execution_time_espt, execution_time_aespt], ["MC", "ESPT", "AESPT"], title="Execution Time", file_name="%s_time"%norad_id)
	Visualize.covariance_ellipsoid([df_mc, df_espt, df_aespt], ["MC", "ESPT", "AESPT"], title="Covariance Ellipsoid", file_name="%s_covariance"%norad_id, type='solid', alpha=0.3)
	Visualize.sample_point_visualize([df_mc, df_espt, df_aespt],  ["MC", "ESPT", "AESPT"], title="Position Sample Points", file_name="%s_samples"%norad_id, alpha=0.3)
	Visualize.comparision_table([df_mc, df_espt, df_aespt],  ["MC", "ESPT", "AESPT"], title1="Comparision between different techniques", title2="%Comparision between different techniques", file_name="%s_table"%norad_id)


	# Visualize.PCA(df_mc.loc[df_mc["timestep"] == 0, ["x", "y", "z"]], "pca_mc")
	# Visualize.PCA(df_espt.loc[df_espt["timestep"] == 0, ["x", "y", "z"]], "pca_espt")
	# Visualize.PCA(df_aespt.loc[df_aespt["timestep"] == 0, ["x", "y", "z"]], "pca_aestp")

	# Visualize.covariance_and_samples(df_mc, file_name="%s_covariance_sample_mc"%norad_id)
	# Visualize.covariance_and_samples(df_espt, file_name="%s_covariance_sample_espt"%norad_id)
	# Visualize.covariance_and_samples(df_aespt, file_name="%s_covariance_sample_aespt"%norad_id)