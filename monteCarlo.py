import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# function to calculate collision probability using monte carlo approach
# states_1_ca and states_2_ca are the orbital states at closest approach                                           collision thershold in meters
def monte_carlo_collision_probability(states_1_ca, cov_matrix_1_ca, states_2_ca, cov_matrix_2_ca, num_samples=100000, collision_threshold=30*u.m):
   
    # Sample from the error distributions
    error_samples_1 = np.random.multivariate_normal([0, 0, 0], cov_matrix_1_ca, size=num_samples)*u.m
    error_samples_2 = np.random.multivariate_normal([0, 0, 0], cov_matrix_2_ca, size=num_samples) *u.m
    # print("error samples")
    # for i in error_samples_1:

    #     print(i)
    # Calculate the relative positions for each sample
    positions_1_samples = states_1_ca.r + error_samples_1
    positions_2_samples = states_2_ca.r + error_samples_2
    
    
    # Calculate the distances between the objects for each sample
    distances_samples = np.linalg.norm(positions_1_samples - positions_2_samples, axis=1)

    

    # Create an array of indices for the x-axis
    indices = np.arange(len(distances_samples))

    # # Plot distance_samples against indices
    # plt.plot(indices, distances_samples, label='Distance Samples')

    # # Plot a straight line at y=50
    # plt.axhline(y=50, color='red', linestyle='--', label='y=50')

    # # Add labels and title
    # plt.xlabel('Index')
    # plt.ylabel('Distance Samples')
    # plt.title('Distance Samples and y=50 Line')

    # # Add a legend
    # plt.legend()

    # # Display the plot
    # plt.show()

    # Count the number of samples where the distance is below the threshold less than 1m
    collision_samples = np.sum(distances_samples < collision_threshold)

    # Calculate the probability of collision
    collision_probability = collision_samples / num_samples

    return collision_probability





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
    s1=5 # radius of satellite

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
    # Example usage:
    collision_prob = monte_carlo_collision_probability(propagated_object1, cov_matrix_1, propagated_object2, cov_matrix_2)
    print(f"Monte Carlo Probability of Collision: {collision_prob:.4f}")
