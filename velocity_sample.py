import numpy as np
from numpy.linalg import norm
from astropy import units as u
import pdb
from astropy.constants import G, M_earth
from poliastro.constants import R_earth
mu = G*M_earth

def theta_rv(nu,e):
    theta = np.arcsin((1 + e*np.cos(nu))/np.sqrt(1 + 2*e*np.cos(nu) + e**2))
    return theta

def vmag(r, nu, e):
    vmag = np.sqrt(mu/r*(2 - (1 - e**2)/(1 + e*np.cos(nu))))
    return vmag

def velocity_unit(r_z, theta):
    r_u = r_z/norm(r_z)
    #const = np.sqrt((r_u[0] - r_u[1])**2 + (r_u[0] + r_u[1])**2 + r_u[2]**2)
    min_cos_inc = -np.sqrt((r_u[0]**2 + r_u[1]**2)/np.sin(theta)**2*(1 - np.cos(theta)**2))
    max_cos_inc = np.sqrt((r_u[0]**2 + r_u[1]**2)/np.sin(theta)**2*(1 - np.cos(theta)**2))
    cos_inc = np.random.uniform(min_cos_inc, max_cos_inc)
    c1 = r_u[1]/r_u[0]
    c2 = np.sin(theta)*cos_inc/r_u[0]
    c3 = -(r_u[0]**2 + r_u[1]**2)/(r_u[0]*r_u[2])
    c4 = -r_u[1]/(r_u[0]*r_u[2])*np.sin(theta)*cos_inc + np.cos(theta)/r_u[2]
    v_roots = np.roots([(1 + c1**2 + c3**2), 2*(c1*c2 + c3*c4), (c2**2 + c4**2 -1)])
    V1 = np.array([v_roots[0], c1*v_roots[0] + c2, c3*v_roots[0] + c4])
    V2 = np.array([v_roots[1], c1*v_roots[1] + c2, c3*v_roots[1] + c4])
    return np.array([V1, V2])

def velocity_sample(r_z, min_alt, Ns):
    """_summary_
    Genrates debsris velocity sample for a given satellite position for collision

    Args:
        r_z (numpy array): Satellite position
        min_alt (scaler in meter): minimum allowable altitude
        Ns (scaler): Number of samples to be generated
    """
    
    r = norm(r_z)*1000*u.m
    min_r = R_earth + min_alt*u.m
    print('r', r)
    print('min_r', min_r)
    #pdb.set_trace()
    min_ta = np.arccos((2*min_r - r)/r)
    print('min_ta', min_ta)
    
    true_anomaly = np.random.uniform(min_ta.value, 2*np.pi - min_ta.value, (Ns))*u.rad
    theta = np.zeros((Ns))
    ecc_sample = np.zeros((Ns))
    for (ta, i) in zip(true_anomaly, range(0,Ns)):
        e_max = (r  - min_r)/(min_r - r*np.cos(ta))
        ecc_sample[i] = np.random.uniform(0,e_max.value)
        theta[i] = theta_rv(ta.value, ecc_sample[i])
    v_unit = np.zeros((2*Ns,3))
    for (th, i) in zip(theta, range(0,Ns)):
        v_unit[i,:] = velocity_unit(r_z, th)[0,:]
        v_unit[Ns + i,:] = velocity_unit(r_z, th)[1,:]
    
    v_mag = vmag(r, np.concatenate((true_anomaly, true_anomaly)), np.concatenate((ecc_sample, ecc_sample)))/1000
    v_sample = np.array([v_unit[i,:]*v_mag[i] for i in range(0, len(v_unit))])
    return v_sample

