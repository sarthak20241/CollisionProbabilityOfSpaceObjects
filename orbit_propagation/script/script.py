import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

import argparse
from dateutil.parser import parse as date_parser

def convertStrToDate(date):
    return date_parser(date).date()

all_techniques = ["MC", "ESPT", "AESPT"]

parser =  argparse.ArgumentParser(prog="Sample propagation", description="Propagates samples using different techniques")


parser.add_argument("-id", "--norad", nargs='+', help="Norad IDs of objects")
parser.add_argument("-ttle", "--tle-time-range", dest="tle_time_range", nargs=2, help="Time range to collect TLE data from(dd/mm/yyyy)", type=convertStrToDate)
parser.add_argument("-n", "--samples", type=int, default=100, help="Specify number of sample points to generate")
parser.add_argument("-tprop", "--propagate-by", nargs='+', dest="propagateby", type=int, help="Timesteps to propagate by in seconds")
parser.add_argument("-t", "--technique", choices=all_techniques, nargs='+', default=all_techniques, help="Technique for propagation")

parser.add_argument("-v", "--verbose", action="store_true", help="Print output at every step")
parser.add_argument("-p", "--plot", action="store_true", help="Plot propagation graphs")
args = parser.parse_args()

from OUP.TLEData import TLEFetch
from OUP import ParticleSampling
from OUP.TLEDifferentiate import TLEDifferentiate
from OUP.ParticleSampling import ParticleSampling
from OUP.ESPT import ESPT
from OUP.MC import MC
from OUP.Visualize import Visualize

#TLE Fetching
fetcher_obj = TLEFetch("sanat@iiitd.ac.in", "DhgiQQvdVlOfHT49", verbose=args.verbose)
data = fetcher_obj.get_data(args.norad, time_range=args.tle_time_range)

#TLE differentiating
diff_obj = TLEDifferentiate(args.verbose)
norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, norad_rv = diff_obj.differentiate(data, save=True)

particle_obj = ParticleSampling(True)
norad_samples = particle_obj.generate_samples(norad_mean_r, norad_mean_v, norad_covariance_r, norad_covariance_v, data.get_norads(), args.samples)


mc_obj = MC(args.verbose, False)
espt_obj = ESPT(args.verbose, False)



if args.plot:
    for (norad_id, sample_points) in norad_samples.items():
        li = [None]*len(args.technique)
        if "MC" in args.technique:
            mc = mc_obj.propagate(sample_points, args.propagateby, norad_ids=norad_id)
            li[args.technique.index("MC")] = mc
            
        if "ESPT" in args.technique:
            espt = espt_obj.propagate(sample_points, args.propagateby, norad_id=norad_id)
            li[args.technique.index("ESPT")] = espt
            
        if "AESPT" in args.technique:
            aespt = espt_obj.propagate(sample_points, args.propagateby, norad_id=norad_id, propagation_type="aespt")
            li[args.technique.index("AESPT")] = aespt

        
        Visualize.execution_time([info[2] for info in li], args.technique, file_name="%s_time"%norad_id)
        Visualize.covariance_ellipsoid([info[0] for info in li], args.technique, title="%s_covariance"%norad_id, file_name="%s_covariance"%norad_id)
        Visualize.sample_point_visualize([info[0] for info in li],  args.technique, file_name="%s_samples"%norad_id)
        
        for i in range(len(li)):
            Visualize.covariance_and_samples(li[i][0], file_name="%s_covariance_sample_%s"%(norad_id, args.technique[i]))