import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
import numpy as np
from scipy.optimize import minimize

def projectOnEncounterPlane(coords,U):
    return U@coords

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




# function to calculate tranformation matrix
def calculate_U(v_s, v_d):
    # Calculate the relative velocity vector
    v_r = v_d - v_s
    
    # Normalize the vectors
    #basis vectors of encounter plane
    i_hat = v_r / np.linalg.norm(v_r)
    j_hat = np.cross(v_d, v_s) / np.linalg.norm(np.cross(v_d, v_s))
    k_hat = np.cross(i_hat, j_hat)
    
    # Construct the transformation matrix C
    U = np.column_stack((i_hat, j_hat, k_hat))
    Uinv = np.linalg.inv(U)
    
    return Uinv



def calculate_T(v_r):
    # Calculate the new unit vectors using the given formula


    i_new = v_r / np.linalg.norm(v_r)
    j_new = np.cross(v_d, v_s) / np.linalg.norm(np.cross(v_d, v_s))
    k_new = np.cross(i_new, j_new)

    # Create the transformation matrix T    
    # basis vectors for encounter plane
    T = np.vstack([i_new, j_new, k_new])

    return T

def project_to_conjunction_plane(P_3d, T):
    # Project the 3D error covariance matrix into the conjunction plane
    P_2d_conjunction = T @ P_3d @ T.T

    return P_2d_conjunction




# function to calculate the collision parametes that would be require to calculate collision probability 
# U- Tranformation Matrix
#sigma- standard deviation
#q - miss vector - relative displacement between both objects
def calculate_collision_parameters(U, sigma_x, sigma_y, sigma_z, q):
    """
    Returns:
    - a, b, c, d, e, f, g, alpha, beta: Parameters from the research paper
    - T: Rotation matrix
    - qr: Transformed relative displacement vector
    """
    # Calculate elements of the transformation matrix U
    U11, U12, U13 = U[0]
    U21, U22, U23 = U[1]
    U31, U32, U33 = U[2]

    # Calculate parameters a, b, c, d, e, f, g
    a = (U13**2) / (2 * sigma_x**2) + (U23**2) / (2 * sigma_y**2) + (U33**2) / (2 * sigma_z**2)
    c = (U11 * U13) / sigma_x**2 + (U21 * U23) / sigma_y**2 + (U31 * U33) / sigma_z**2
    d = (U12 * U13) / sigma_x**2 + (U22 * U23) / sigma_y**2 + (U32 * U33) / sigma_z**2
    e = (U11**2) / (2 * sigma_x**2) + (U21**2) / (2 * sigma_y**2) + (U31**2) / (2 * sigma_z**2) - (c**2) / (4 * a)
    f = (U12**2) / (2 * sigma_x**2) + (U22**2) / (2 * sigma_y**2) + (U32**2) / (2 * sigma_z**2) - (d**2) / (4 * a)
    g = (U11 * U12) / sigma_x**2 + (U21 * U22) / sigma_y**2 + (U31 * U32) / sigma_z**2 - (c * d) / (2 * a)
    

    alpha = (e + f) / 2 - np.sqrt(g**2 + (f - e)**2) / 2
    beta = (e + f) / 2 + np.sqrt(g**2 + (f - e)**2) / 2

    # Calculate rotation angle phi
    cos_phi = np.sqrt(0.5 * (1 - (f - e) / np.sqrt(g**2 + (f - e)**2)))
    sin_phi_positive = np.sqrt(0.5 * (1 + (f - e) / np.sqrt(g**2 + (f - e)**2)))
    sin_phi_negative = -sin_phi_positive


    # Choose the sign to make the cross term zero
    cross_term = 2 * (f - e) * sin_phi_negative * cos_phi + g * (cos_phi**2 - sin_phi_negative**2)

    if np.isclose(cross_term, 0):
        sin_phi = sin_phi_negative
    else:
        sin_phi = sin_phi_positive

    
    # Calculate rotation matrix T
    T = np.array([[cos_phi, sin_phi],
                  [-sin_phi, cos_phi],])

    

    # Calculate transformed relative displacement vector qr
    qr = np.dot(T, q)

    return a, c, d, e, f, g, alpha, beta, T, qr


#scales the y parameter of qr to produce qrs
def qrscaleY(qr, alpha, beta):
   
    # Calculate scale factor
    scale_factor = np.sqrt(beta / alpha)

    # Update qrs_2 value
    qrs_2_updated = np.array([0.0, scale_factor * qr[1]])

    return qrs_2_updated




def rotate_vector(vector, angle):
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle), np.cos(angle)]])
    return np.dot(rotation_matrix, vector)

#function to calculate collision probability using patera method
# alpha beta are the parameters derived using Transformation matrix
# q_r and q_r_s are transformed and scaled miss vector respectively
# s is the combined radius of both objects /hard body circle
def calculate_collision_probability(alpha, beta, q_r, q_r_s, s,rprimary,rsecondary, steps=1000):
    prob_sum = 0
    theta = 2 * np.pi  / steps
    X = np.array([1, 0])
    for i in range(steps):
        
        
        X_prime = rotate_vector(X, theta)
        
        X_m_prime = np.dot(np.array([[s, 0], [0, s * np.sqrt(beta / alpha)]]), X_prime)
        X_e_prime = X_m_prime + np.array([q_r[0], q_r_s[1]])

        X_e = np.dot(np.array([[s, 0], [0, s * np.sqrt(beta / alpha)]]), X) + np.array([q_r[0], q_r_s[1]])

        d_theta = np.arcsin(np.cross(X_e, X_e_prime) / (np.linalg.norm(X_e) * np.linalg.norm(X_e_prime)))

        int_value = np.exp(-alpha * np.dot((X_e + X_e_prime) / 2, (X_e + X_e_prime) / 2))
        prob_sum += int_value * d_theta

        X=X_prime

    dist = np.linalg.norm(rsecondary)
    if(dist<s):
        prob=(1 - prob_sum)/(2* np.pi)
    else:
        prob = -prob_sum / (2 * np.pi)
    
    if prob>1:
        prob= 1-prob
    if prob<0:
        prob=1+prob

    # prob = -prob_sum / (2 * np.pi)

    return prob






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

    v_s=propagated_object1.v/(u.km/u.s)
    v_d=propagated_object2.v/(u.km/u.s)
    # v_s = np.array([1.0, 0.0, 0.0])  # Example relative velocity vector of the space station
    # v_d = np.array([0.0, 1.0, 0.0])  # Example relative velocity vector of the debris

    # transformation matrix from 3-d plane to encounter plane
    U = calculate_U(v_s, v_d)



    # for short term encounters relative velocity constant => error covariance matrix could be added to get relative error
    error_covariance_relative = cov_matrix_1+cov_matrix_2



    sigma_x, sigma_y, sigma_z = get_standard_deviations(error_covariance_relative)

    s=s1+s2 # radius of initial hard body circle
    #q = ((propagated_object1.r-propagated_object2.r)/u.km)*1000
    q = np.array([100, 0]) # second object's location (hard body circle) relative displacement vector 
    # hard body circle centered at q
    # the origin might be on the center of covariance ellipse of primary object

    # Call the function
    a, c, d, e, f, g, alpha, beta, T, qr = calculate_collision_parameters(U, sigma_x, sigma_y, sigma_z, q)



    # Call the function
    qrs = qrscaleY(qr, alpha, beta)




    collision_probability = calculate_collision_probability(alpha, beta, qr, qrs, s,r0_object1,r0_object2)
    print("Patera Collision Probability:", collision_probability)