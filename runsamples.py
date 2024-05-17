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
        


# Initial position and velocity vectors
r0_object1 = np.array([7500, 500, 320]) # in km
v0_object1 = np.array([12, 2, 2]) # in km/s

v0_obj1=getPrimaryVelocity(r0_object1)


# vsamp=vs.velocity_sample(r0_object1,500,100)
# print(vsamp)

monte_carlo=np.array([])
alfano=np.array([])
patera=np.array([])
deltaR=np.array([])

for i in range(15):
    print("\n\n Sample ",i+1)
    print("\n")

    mnc,alf,pat=r.compareResults(0)
    # mnc,alf,pat=r.compareResults(0.01*(i+1))
    monte_carlo=np.append(monte_carlo,mnc)
    alfano=np.append(alfano,alf)
    patera=np.append(patera,pat)
    
    deltaR=np.append(deltaR,0.01*i)

# Plotting the data
plt.plot(deltaR, monte_carlo, label='Monte Carlo')
plt.plot(deltaR, alfano, label='Alfano')
plt.plot(deltaR, patera, label='Patera')

# Adding labels and title
plt.xlabel('DeltaR(in km)')
plt.ylabel('Collision Probability')
plt.title('Collision Probability vs DeltaR')

# Adding a legend
plt.legend()

# Save the plot as a vector image ( SVG)
plt.savefig('collision_plot.svg', format='svg') 
# Display the plot
plt.show()





    