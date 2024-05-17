from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
import numpy as np
from scipy.optimize import minimize


# Function to calculate TCPA using the formula from the research paper
def calculate_tcpa(r0_s, v0_s, r0_d, v0_d):
    # Calculate relative position vector and relative velocity vector
    rho_o = r0_s - r0_d
    v_r = v0_s - v0_d

    print("rho_0 ",rho_o)
    print("v_r",v_r)
    # Calculate tcpa
    tcpa = -rho_o.dot(v_r) / v_r.dot(v_r)
    return tcpa


if __name__ == "__main__":



    # Initial position and velocity vectors
    r0_object1 = [7000, 0, 0] * u.km
    v0_object1 = [76, 7, 0] * u.km / u.s
    s1=5 # radius of satellite 100kg class microsattelite

    r0_object2 = [7100, 0, 0] * u.km
    v0_object2 = [65, 7.1, 0] * u.km / u.s
    s2=0.7 # radius of debris


    # Define initial orbits
    orbit_object1 = Orbit.from_vectors(Earth, r0_object1, v0_object1)
    orbit_object2 = Orbit.from_vectors(Earth, r0_object2, v0_object2)




    # Calculate TCPA
    tcpa = calculate_tcpa(r0_object1, v0_object1, r0_object2, v0_object2)
    print("Time of Closest Approach (TCPA):", tcpa)

    # Propagate orbits to TCPA
    propagated_object1 = orbit_object1.propagate(tcpa)
    propagated_object2 = orbit_object2.propagate(tcpa)

    cov_matrix_1=np.diag([13, 15, 9])
    cov_matrix_2=np.diag([16, 23, 8])

    # Print propagated orbits
    print("\nPropagated Orbit for Object 1 at TCPA:")
    print(propagated_object1.r)
    print("\nPropagated Orbit for Object 2 at TCPA:")
    print(propagated_object2.r)

    miss_vector=propagated_object1.r-propagated_object2.r
    print("Miss Vector ",miss_vector)
    v_r=propagated_object1.v=propagated_object2.v
    print("Relative Velocity",v_r)
