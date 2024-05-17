import alfano as af
import monteCarlo as mc
import patera as pt
import calculateTcpa as tc
import run as r

from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

import velocity_sample as vs

def getPrimaryVelocity(r):
    vsamp=vs.velocity_sample(r,r[1]-100,1)
    return vsamp[0]
def getSecondaryVelocity(r,vel):
    v=None
    
    # vsamp=vs.velocity_sample(coords,500,100)
    # print(vsamp)
    while v is None:
        
        vsamp=vs.velocity_sample(r,r[1]-100,100)
        
        for i in vsamp:
            
            print(i)
            x=i[0]-vel[0]
            y=i[1]-vel[1]
            z=i[2]-vel[2]
            if(x**2+y**2+z**2>5):
                v=i
                break
    
    return v
        
def generate_random_r1_LEO():
    # Random altitude within LEO range (in kilometers)
    altitude_km = np.random.uniform(160, 2000)
    
    
    
    # Generate random angle theta (longitude) in radians
    theta = np.random.uniform(0, 2*np.pi)
    
    # Generate random position vector in the equatorial plane
    r1 = [altitude_km * np.cos(theta), altitude_km * np.sin(theta), 0]
    
    return r1

def sureShotSimulation(r1,v1,sat_name):
    # Initial position and velocity vectors
    r0_object1 = r1 / (np.array([1,1,1]) *u.km)# in km
    v1 = v1 / (np.array([1,1,1])*u.km/u.s)# in km/s


    # vsamp=vs.velocity_sample(r0_object1,500,100)
    # print(vsamp)

    monte_carlo=np.array([])
    alfano=np.array([])
    patera=np.array([])
    sigma=np.array([])

    for sig in range(20, 100):
        print("\n\n Sample ",sig-19)
        print("\n")
        
        r1=np.array([7500,1500,200])
        print("Position vector of object 1 ", r1)
        mnc,alf,pat=r.compareSureShotResults2(sig,r1,v1)
        
        monte_carlo=np.append(monte_carlo,mnc)
        alfano=np.append(alfano,alf)
        patera=np.append(patera,pat)
        
        sigma=np.append(sigma,sig)

    # Plotting the data
    plt.plot(sigma, monte_carlo, label='Monte Carlo')
    plt.plot(sigma, alfano, label='Alfano')
    plt.plot(sigma, patera, label='Patera')

    # Adding labels and title
    plt.xlabel('Standard Deviation in x,y,z direction (in m)')
    plt.ylabel('Collision Probability')
    plt.title('Collision Probability vs Standard Deviation for '+sat_name)

    # Adding a legend
    plt.legend()

    # Save the plot as a vector image ( SVG)
    plt.savefig('sure_shot_collision_plot'+sat_name+'.png', format='png') 
    # Display the plot
    plt.close()
    print()

if __name__ == "__main__":
    # Initial position and velocity vectors
    r0_object1 = np.array([7500, 500, 320]) # in km
    v0_object1 = np.array([12, 2, 2]) # in km/s

    v1=getPrimaryVelocity(r0_object1)


    # vsamp=vs.velocity_sample(r0_object1,500,100)
    # print(vsamp)

    monte_carlo=np.array([])
    alfano=np.array([])
    patera=np.array([])
    sigma=np.array([])

    for sig in range(20, 501):
        print("\n\n Sample ",sig-19)
        print("\n")
        
        r1=np.array([7500,1500,200])
        print("Position vector of object 1 ", r1)
        mnc,alf,pat=r.compareSureShotResults2(sig,r1,v1)
        
        monte_carlo=np.append(monte_carlo,mnc)
        alfano=np.append(alfano,alf)
        patera=np.append(patera,pat)
        
        sigma=np.append(sigma,sig)

    # Plotting the data
    plt.plot(sigma, monte_carlo, label='Monte Carlo')
    plt.plot(sigma, alfano, label='Alfano')
    plt.plot(sigma, patera, label='Patera')

    # Adding labels and title
    plt.xlabel('Sigma(in m)')
    plt.ylabel('Collision Probability')
    plt.title('Collision Probability vs Sigma')

    # Adding a legend
    plt.legend()

    # Save the plot as a vector image ( SVG)
    plt.savefig('sure_shot_collision_plot.svg', format='svg') 
    # Display the plot
    plt.show()