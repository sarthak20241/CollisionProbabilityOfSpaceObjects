from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
from poliastro.util import norm
import numpy as np
import missCollision2 as ms
import sureShot2 as ss
import OUP.TLEData as tle
import datetime


# Define an empty dictionary to store the orbits
orbits = {}
missDebris = {}

time_range_tle = (datetime.date(2023, 4, 30), datetime.date(2023, 5, 1))
norad_ids = ["42962", "46826", "43286", "45026"]

#TLE Fetching
fetcher_obj = tle.TLEFetch("sanat@iiitd.ac.in", "DhgiQQvdVlOfHT49")
data = fetcher_obj.get_data(norad_ids, time_range=time_range_tle)


TLE_data = [("Iridium 136",data['42962'][0][0],data['42962'][0][1]),("NAVSTAR 80",data['46826'][0][0],data['46826'][0][1]),("GSAT-30",data['45026'][0][0],data['45026'][0][1]),("IRNSS-1I",data['43286'][0][0],data['43286'][0][1])]
print(TLE_data)
def main():
    # Parse TLE and define orbits for all satellites
    for name, line1, line2 in TLE_data:
        satellite = twoline2rv(line1, line2, wgs84)
        
        position, velocity = satellite.propagate(2024, 3, 25, 20, 0, 0)  # Note: Time is in UTC
        
        position = position * u.km
        velocity = velocity * (u.km / u.s)
        # print(velocity)
        orbit = Orbit.from_vectors(Earth, position, velocity)
        
        orbits[name] = orbit

    # Print all orbits
    for name, orbit in orbits.items():
        print(f"{name} orbit:")
        print(orbit)
        print()

        debris=ms.missCollisionSimulations(orbit.r,orbit.v,name)
        missDebris[name]=debris
        ss.sureShotSimulation(orbit.r,orbit.v,name)

if __name__ == '__main__':

    main()
