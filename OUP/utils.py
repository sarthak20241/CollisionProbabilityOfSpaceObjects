import datetime
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import astropy.units as u 

def date_range_list(time_range, date_format="%Y-%m-%d"):
        (start, end) = time_range
        return [(start+datetime.timedelta(days=day)).strftime(date_format) for day in range((end-start).days+1)]

def create_poliastro_object(r, v):
        r = r.values.tolist()*u.km
        v = v.values.tolist()*u.km/u.s
        return Orbit.from_vectors(Earth,r, v)