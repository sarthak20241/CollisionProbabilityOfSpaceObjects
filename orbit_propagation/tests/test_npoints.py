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
norad_objs = {"42962":{"code":"LEO", "timesteps":range(0, 6001, 600)}, "46826":{"code":"MEO", "timesteps":range(0, 43081, 4308)}, "43286":{"code":"GSO", "timesteps":range(0, 86161, 8616)}, "45026":{"code":"GEO", "timesteps":range(0, 86161, 8616)}}
num_sample_points = [100, 500, 1000]

#TLE Fetching
fetcher_obj = TLEFetch("sanat@iiitd.ac.in", "DhgiQQvdVlOfHT49")
data = fetcher_obj.get_data(list(norad_objs.keys()), time_range=time_range_tle)

#TLE differentiating
diff_obj = TLEDifferentiate()
norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_rv = diff_obj.differentiate(data, save=True)

particle_obj = ParticleSampling(True)
mc_obj = MC(False, False)
espt_obj = ESPT(verbose=False, save=False)

for n in num_sample_points:	
    norad_samples = particle_obj.generate_samples(norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, data.get_norads(), n)
    for (norad_id, sample_points) in norad_samples.items():
        df_mc, mean_df_mc, execution_time_mc = mc_obj.propagate(sample_points, norad_objs[norad_id]["timesteps"], norad_ids=norad_id)
        df_espt, mean_df_espt, execution_time_espt = espt_obj.propagate(sample_points, norad_objs[norad_id]["timesteps"], norad_id=norad_id)
        df_aespt, mean_df_aespt, execution_time_aespt = espt_obj.propagate(sample_points, norad_objs[norad_id]["timesteps"], norad_id=norad_id, propagation_type="aespt")

        code = norad_objs[norad_id]["code"]
        
        Visualize.execution_time([execution_time_mc, execution_time_espt, execution_time_aespt], ["MC", "ESPT", "AESPT"], title="Execution time for %s with %s sample points"%(code, n), file_name="%s_%s_time"%(code, n))
        Visualize.covariance_ellipsoid([df_mc, df_espt, df_aespt], ["MC", "ESPT", "AESPT"], title="%s Covariance for %s sample points"%(code, n), file_name="%s_%s_covariance"%(code, n), type='solid', alpha=0.3)
        Visualize.sample_point_visualize([df_mc, df_espt, df_aespt],  ["MC", "ESPT", "AESPT"], title="Sample points for %s with %s sample points"%(code, n), file_name="%s_%s_samples"%(code, n), alpha=0.3)
        Visualize.comparision_table([df_mc, df_espt, df_aespt],  ["MC", "ESPT", "AESPT"], title1="%s: %s points"%(code, n), title2="%s: %s points"%(code, n), file_name="%s_%s_table"%(code, n))


	# Visualize.PCA(df_mc.loc[df_mc["timestep"] == 0, ["x", "y", "z"]], "pca_mc")
	# Visualize.PCA(df_espt.loc[df_espt["timestep"] == 0, ["x", "y", "z"]], "pca_espt")
	# Visualize.PCA(df_aespt.loc[df_aespt["timestep"] == 0, ["x", "y", "z"]], "pca_aestp")