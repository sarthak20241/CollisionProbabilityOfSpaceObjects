import numpy as np
from scipy.special import erf
from scipy.integrate import simps
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
import numpy as np
from scipy.optimize import minimize
import math

#function to calculate collision probability using alfano's method
#sigma- standard deviation
#xm,ym- x and y component of miss vector
# obj_radius - radius of hard body circle / radius of both objects combined
def collision_probability(sigma_x, sigma_y, xm, ym, obj_radius):

    OBJ = obj_radius
    
    try:
        m = int(5 * OBJ / min(sigma_x, sigma_y, np.sqrt(xm**2 + ym**2)))
        m = max(10, min(50, m))
    except:
        m=50
   

    m_odd = 0
    dx = OBJ / (2 * m)
    for i in range(1, m):
        x = OBJ * (2 * i - m) / m
        y = np.sqrt(OBJ**2 - x**2)

        term1 = erf((ym + y) / (sigma_y * np.sqrt(2)))
        term2 = erf((-ym + y) / (sigma_y * np.sqrt(2)))

        exponential_term1 = np.exp(-((xm + x)**2) / (2 * sigma_x**2))
        exponential_term2 = np.exp(-((xm - x)**2) / (2 * sigma_x**2))

        m_odd += 4 * (term1 + term2) * (exponential_term1+exponential_term2)

    m_0 = 0
    
    x = (0.015 * dx) - OBJ
    y = np.sqrt(OBJ**2 - x**2)

    term1 = erf((ym + y) / (sigma_y * np.sqrt(2)))
    term2 = erf((-ym + y) / (sigma_y * np.sqrt(2)))

    exponential_term1 = np.exp(-((xm + x)**2) / (2 * sigma_x**2))
    exponential_term2 = np.exp(-((xm - x)**2) / (2 * sigma_x**2))

    m_0 += 2 * (term1 + term2) * (exponential_term1+exponential_term2)

    


    m_even=0
    dx = OBJ / (2 * m)
    for i in range(1, m):
        x = ((2 * i) *dx ) - OBJ
        y = np.sqrt(OBJ**2 - x**2)

        term1 = erf((ym + y) / (sigma_y * np.sqrt(2)))
        term2 = erf((-ym + y) / (sigma_y * np.sqrt(2)))

        exponential_term1 = np.exp(-((xm + x)**2) / (2 * sigma_x**2))
        exponential_term2 = np.exp(-((xm - x)**2) / (2 * sigma_x**2))

        term3 = erf((ym + OBJ)/ (sigma_y * np.sqrt(2)))
        term4 = erf((-ym + OBJ)/ (sigma_y * np.sqrt(2)))
        
        exponential_term3 = np.exp(-((xm)**2) / (2 * sigma_x**2))

        m_even += 4 * (term1 + term2) * (exponential_term1+exponential_term2)
        m_even += 2 * (term3 + term4) * (exponential_term3)


    
    probability = (dx / (3*(np.sqrt(8 * np.pi) * sigma_x ))) * (m_0 + m_even + m_odd)

    return probability


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
def get_standard_deviations(error_covariance_matrix):
   
    # Extract diagonal elements
    sigma_x_squared = error_covariance_matrix[0, 0]
    sigma_y_squared = error_covariance_matrix[1, 1]
    sigma_z_squared = error_covariance_matrix[2, 2]

    # Calculate standard deviations
    sigma_x = np.sqrt(sigma_x_squared)
    sigma_y = np.sqrt(sigma_y_squared)
    sigma_z = np.sqrt(sigma_z_squared)

    return sigma_x, sigma_y, sigma_z

# Initial position and velocity vectors
r0_object1 = [7000, 0, 0] * u.km
v0_object1 = [76, 7, 0] * u.km / u.s
s1=5 # radius of satellite

r0_object2 = [7100, 0, 0] * u.km
v0_object2 = [65, 7.1, 0] * u.km / u.s
s2=0.7 # radius of debris


# Define initial orbits
orbit_object1 = Orbit.from_vectors(Earth, r0_object1, v0_object1)
orbit_object2 = Orbit.from_vectors(Earth, r0_object2, v0_object2)



if __name__ == "__main__":
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

    error_covariance_relative = cov_matrix_1+cov_matrix_2

    sigma_x, sigma_y, sigma_z = get_standard_deviations(error_covariance_relative) 

    print(sigma_x,sigma_y)


    # sigma_x and sigma_y length of major and minor axis taken 


    # position of secondary object relative to primary on covariance frame
    xm = abs(miss_vector[0]*1000/u.km)    # x component or horizontal distance
    ym = abs(miss_vector[1]*1000 /u.km)   # y component or vertical distance

    print(xm,ym)
    obj_radius = s1+s2 # radius of combined object
    print(obj_radius)
    collision_prob = collision_probability(sigma_x, sigma_y, xm, ym, obj_radius)
    print(f"Alfano Collision Probability: {collision_prob}")
