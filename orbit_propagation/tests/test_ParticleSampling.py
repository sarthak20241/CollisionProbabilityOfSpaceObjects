import sys
import os

# Add the parent directory to the sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from OUP.ParticleSampling import ParticleSampling
from OUP.TLEData import TLEFetch
from OUP.TLEDifferentiate import TLEDifferentiate

norad_ids = ["42962", "46826", "43286", "45026"]
fetcher = TLEFetch("sanat@iiitd.ac.in", "DhgiQQvdVlOfHT49")
data = fetcher.get_data(norad_ids)
diff = TLEDifferentiate()
norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_rv = diff.differentiate(data, save=True)

particle_obj = ParticleSampling(True)
sample_dict = particle_obj.generate_samples(norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_ids, 100)
print(sample_dict)
print(sample_dict[norad_ids[0]])