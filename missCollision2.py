import numpy as np
import run as r
import matplotlib.pyplot as plt
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.util import norm
from astropy import units as u
import velocity_sample as vs
import perpendicularVec as pv






def generate_random_r1_LEO():
    # Random altitude within LEO range (in kilometers)
    altitude_km = np.random.uniform(160, 2000)
    
    # Generate random angle theta (longitude) in radians
    theta = np.random.uniform(0, 2*np.pi)
    
    # Generate random position vector in the equatorial plane
    r1 = [altitude_km * np.cos(theta), altitude_km * np.sin(theta), 0]
    
    return r1


def calculate_v1(r1, v):
    # Ensure r1 is a numpy array
    r1 = np.array(r1)
    
    # Check if r1 is in the equatorial plane
    if r1[2] != 0:
        raise ValueError("Position vector r1 is not in the equatorial plane.")
    
    # Calculate v1 with the desired magnitude v
    v1 = np.array([-r1[1], r1[0], 0])
    
    # Normalize v1
    v1_normalized = v1 / np.linalg.norm(v1)
    
    # Adjust magnitude of v1
    v1 = v * v1_normalized
    
    return v1

def calculate_velocity_magnitude(r):
    """
    Calculate the magnitude of the velocity of an object in a circular orbit.
    """
    # Define constants
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    M = 5.972e24     # Mass of the Earth in kg (example)
    
    # Convert position vector from km to meters
    r_m = r * 1000  # Convert km to meters
    
    # Calculate the magnitude of the position vector in meters
    radius = np.linalg.norm(r_m)
    
    # Calculate the magnitude of the acceleration due to gravity at the orbital radius
    acceleration_due_to_gravity = G * M / radius**2
    
    # Calculate the magnitude of the velocity using the centripetal acceleration formula
    velocity_magnitude = np.sqrt(acceleration_due_to_gravity * radius) / 1000  # Convert m/s to km/s
    
    return velocity_magnitude


def generate_r2(r1, d):
    # Ensure r1 is a numpy array
    r1 = np.array(r1)

    # Generate random spherical coordinates
    theta = np.random.uniform(0, 2 * np.pi)  # Azimuthal angle
    phi = np.arccos(np.random.uniform(-1, 1))  # Polar angle

    # Convert spherical coordinates to Cartesian coordinates
    x = r1[0] + d * np.sin(phi) * np.cos(theta)
    y = r1[1] + d * np.sin(phi) * np.sin(theta)
    z = r1[2] + d * np.cos(phi)

    # Define the position vector of the secondary object
    r2 = np.array([x, y, z])

    
    return r2

def calculate_r2(r1, d):
    # Ensure r1 is a numpy array
    r1 = np.array(r1)
    
    # Calculate the magnitude of r1
    magnitude_r1 = np.linalg.norm(r1)
    
    # Calculate the unit vector in the direction of r1
    if magnitude_r1 != 0:
        u_r1 = r1 / magnitude_r1
    else:
        raise ValueError("Position vector r1 cannot be a zero vector.")
    
    # Calculate r2
    r2 = r1 + d * u_r1
    
    return r2

def generate_v2(r1,r2,v1):

    v2=None
    
    
    while v2 is None:
        
        possible_velocity_vectors=vs.velocity_sample(r2,r2[1]-100,3000)

            
        # Choose a velocity vector from the list that satisfies the condition
        for v in possible_velocity_vectors:
            # Calculate the angle between the vectors
            cos_angle = np.dot(r2 - r1, v - v1) / (norm(r2 - r1) * norm(v - v1))
            angle = np.arccos(cos_angle)
            
            # Convert the angle to degrees
            angle_degrees = np.degrees(angle)/ 1 
            print(angle_degrees)
            if -90*u.deg < angle_degrees < 90 * u.deg:
                v2=v
                break
        
    return np.array(v2)

def isMiss(r1,r2,v1,v):
    # Calculate the angle between the vectors
    cos_angle = np.dot(r2 - r1, v - v1) / (norm(r2 - r1) * norm(v - v1))
    angle = np.arccos(cos_angle)
    
    # Convert the angle to degrees
    angle_degrees = np.degrees(angle)/ 1 
    print(angle_degrees)
    if -90*u.deg < angle_degrees < 90 * u.deg:
        return True
    return False

def generatePerpendicular_v2(r1,r2,v1):
    pvs=pv.perpendicular_vectors(r2,520)
    v2mag=calculate_velocity_magnitude(r2)
    v2res = None
    for v in pvs:
        v2 = v * v2mag
        if(isMiss(r1,r2,v1,v2)):
            v2res=v2
            break
    return v2res
def calculate_v2(r2, v_prime):
    # Ensure r2 is a numpy array
    r2 = np.array(r2)
    
    
    # Calculate v2 with the desired magnitude v_prime
    v2 = np.array([0, 0, v_prime])
    
    return v2


# input miss distance desired magnitude of v1 and v2
def generateSample(d,v1mag,v2mag):
    r1 = generate_random_r1_LEO() # position vector in equitorial plane in leo orbit (in kms)
    print("Random position vector r1 in LEO:", r1)
    r2 = calculate_r2(r1, d)    # position vector of object 2 d distance away (in kms)
    print("Position vector r2:", r2)
    v1 = calculate_v1(r1, v1mag)    # velocity vector of object 1 (in km/s)
    print("Velocity vector v1:", v1)
    v2 = calculate_v2(r2, v2mag) # velcity vector of object 2 (in km/s)
    print("Velocity vector v2:", v2)

    return r1, v1, r2, v2

def generate_secondary_rv(r1,v1,d):

    
    r1wu=r1/(np.array([1,1,1])*u.km) # r1 without units
    v1wu=v1/(np.array([1,1,1])*(u.km/u.s)) # v1 without units
    v2 = None
    while v2 is None:
        # Calculate the position vector of the secondary object
        r2 = generate_r2(r1,d) 
        v2 = generatePerpendicular_v2(r1wu,r2,v1wu)
    v2 = v2 * (u.km/u.s) 
    # v2 = generate_v2(r1wu,r2,v1wu)* (u.km/u.s)
    r2 = r2 *u.km
    return r2,v2

def missCollisionSimulations(r1,v1,sat_name):
    monte_carlo=np.array([])
    alfano=np.array([])
    patera=np.array([])
    distance=np.array([])
    orbits = []
    # Loop over distances from 0.005 km to 350m in 75 steps
    for i in range(75):
        d = 0.005 * (i+1)  # Distance from 0.005 km to 350m in 75 steps
        r2,v2=generate_secondary_rv(r1,v1,d)
        orbit = Orbit.from_vectors(Earth, r2, v2)
        
        # Store the orbit in the list
        orbits.append(orbit)
        s1=10
        s2=1
        #in meter square
        cov1=np.diag([40000, 40000, 40000])
        cov2=np.diag([200**2, 200**2, 200**2])
        print("\n\n Sample ",i+1)
        print("\n")

        mnc,alf,pat=r.compareMissResults(r1,r2,v1,v2,s1,s2,cov1,cov2)
    
        monte_carlo=np.append(monte_carlo,mnc)
        alfano=np.append(alfano,alf)
        patera=np.append(patera,pat)
        distance=np.append(distance,d)
        

    # Plotting the data
    plt.plot(distance, monte_carlo, label='Monte Carlo')
    plt.plot(distance, alfano, label='Alfano')
    plt.plot(distance, patera, label='Patera')

    # Adding labels and title
    plt.xlabel('Miss Distance (in Km)')
    plt.ylabel('Collision Probability')
    plt.title('Collision Probability vs Miss Distance for '+sat_name)

    # Adding a legend
    plt.legend()

    # Save the plot as a vector image ( SVG)
    plt.savefig('miss_collision_plot'+sat_name+'.png', format='png') 
    plt.close()
    # Display the plot
    #plt.show()
    print()
    return orbits

if __name__ == "__main__":

    monte_carlo=np.array([])
    alfano=np.array([])
    patera=np.array([])
    distance=np.array([])

    # Loop over distances from 0.001 km to 1 km in 1000 steps
    for i in range(150):
        d = 0.005 * (i+1)  # Distance from 0.001 km to 1 km in 1000 steps
        v1mag = np.random.uniform(10, 50)  # Random magnitude for v1
        v2mag = np.random.uniform(10, 50)  # Random magnitude for v2
        r1, v1, r2, v2 = generateSample(d, v1mag, v2mag)
        r1 = r1* u.km 
        r2 = r2* u.km 
        v1 = v1 * (u.km/u.s)
        v2 = v2 * (u.km/u.s)
        s1=10
        s2=1
        #in meter square
        cov1=np.diag([200**2, 200**2, 200**2])
        cov2=np.diag([200**2, 200**2, 200**2])
        print("\n\n Sample ",i+1)
        print("\n")

        mnc,alf,pat=r.compareMissResults(r1,r2,v1,v2,s1,s2,cov1,cov2)
    
        monte_carlo=np.append(monte_carlo,mnc)
        alfano=np.append(alfano,alf)
        patera=np.append(patera,pat)
        distance=np.append(distance,d)
        

    # Plotting the data
    plt.plot(distance, monte_carlo, label='Monte Carlo')
    plt.plot(distance, alfano, label='Alfano')
    plt.plot(distance, patera, label='Patera')

    # Adding labels and title
    plt.xlabel('Miss Distance (in m)')
    plt.ylabel('Collision Probability')
    plt.title('Collision Probability vs Miss Distance')

    # Adding a legend
    plt.legend()

    # Save the plot as a vector image ( SVG)
    plt.savefig('miss_collision_plot.svg', format='png') 
    # Display the plot
    plt.show()