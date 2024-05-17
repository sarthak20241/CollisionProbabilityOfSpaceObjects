import sys
import os

# Add the parent directory to the sys.path
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

norad_ids = ["42962", "46826", "43286", "45026"]
fetcher = TLEFetch("sanat@iiitd.ac.in", "DhgiQQvdVlOfHT49")
data = fetcher.get_data(norad_ids)
diff = TLEDifferentiate()
norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_rv = diff.differentiate(data, save=True)

particle_obj = ParticleSampling(True)
sample_dict = particle_obj.generate_samples(norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_ids, 100)

mc_obj = MC(False, True)
mc_prop_dict = mc_obj.propagate(sample_dict, [0, 10], 1)


espt_obj = ESPT(verbose=True, save=True)
for (norad_id, sample_points) in sample_dict.items():
	df, mean_df, execution_time = espt_obj.propagate(sample_points, [0, 10], norad_id=norad_id)
	break

Visualize.execution_time([mc_prop_dict["42962"][2], execution_time], legend=["MC", "ESPT"])