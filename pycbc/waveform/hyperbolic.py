import numpy as np
from pycbc.types.timeseries import TimeSeries
from numpy import arctanh, tan, sin, cos, array
from scipy.optimize import bisect
import warnings
import pycbc.waveform
warnings.filterwarnings("ignore")

# GR units
pc=3.08*10**(16)
G=6.67*10**(-11)
c=3*10**8
Ms=1.98*10**(30)

## get orbital params
def get_orbit(m1, m2, ri, Qi, vi = np.sqrt(G*(5*10**5)*Ms/(3* 10*pc))):
    
    m1 *= Ms
    m2 *= Ms
    ri *= pc
    Qi = 10**Qi
    
    M = m1+m2
    Mu = m1*m2/(m1+m2)
    GM = G*M
    L = ri*vi*np.sin(Qi)

    phi0 = np.arctan( (L*vi*np.cos(Qi))/(L**2/ri - GM) )
    keep = True
    while keep:
        if phi0 < 0:
            phi0 += np.pi
        elif phi0 > np.pi:
            phi0 -= np.pi
        if phi0 >= 0 and phi0 <= np.pi:
            keep = False
            
    rmin = (L**2)*np.cos(phi0)/(L**2/ri + GM*(np.cos(phi0) - 1))
    p = (L**2)/GM
    ecc = (L**2)/(GM*rmin)-1
    
    return Mu, GM, L, phi0, rmin, p, ecc

# get equation for finding phi as function of t (numerically), and vica-verca    
def Equation(phi, phi0, Ltf, C1, C2, ecc):
    
    PHI = phi0 - phi
    return arctanh( C1* tan(PHI/2)) - C2*sin(PHI)/(1+ecc*cos(PHI)) - Ltf


## Return polarizations
def H_pols(phit, phi0, Mu, GM, D, L, pm, rmin):
    
    factor = 2*G*((GM**2)*Mu)/(D*(L**2)*(c**4))
    PHI = phi0-phit

    trigs = 3*sin(phi0 - 3*phit) + sin(phi0+phit)
    part1 = -6*sin(2*phit)*(1+pm*cos(PHI))**2
    part2 = (pm**2)*(3*sin(2*phit))*(sin(PHI))**2
    part3 = -pm*(1+pm*cos(PHI))*trigs/2
    H12 = factor*(part1+part2+part3)
    
    trigs = 9*cos(-phi0 + 3*phit) - cos(phi0)*cos(phit) + 5*sin(phi0)*sin(phit)
    part1 = -6*cos(2*phit)*(1+pm*cos(PHI))**2
    part2 = (pm**2)*(1+3*cos(2*phit))*(sin(PHI))**2
    part3 = pm*(1+pm*cos(PHI))*trigs/2
    H11 = factor*(part1+part2+part3)
    
    trigs = 9*cos(-phi0 + 3*phit) - 5*cos(phi0)*cos(phit) + sin(phi0)*sin(phit)
    part3 = -pm*(1+pm*cos(PHI))*trigs/2
    H22 = factor*(- part1 - part2+part3)
        
    return array(H11), array(H22), array(H12)
    
def hyperbolic_td(**args):
    
    """
    Docstring:
    Implementation of Sajal's MATHEMATICA code for waveform generation on Python
    
    m1 (type: float):                                          Primary Mass in solar mass units
    m2 (type: float):                                          Secondary Mass in solar mass units
    ri (type: float):                                          initial distance of approach
    Qi (type: float):                                          initial angle of approach in log10 value
    Dl (type float, Default = 6):                              Luminoscity deiscae (in log10 pc values)
    duration (type: float, default = 20):                      time duration
    delta_t (type: float, default = 1/256):                    Sampling rate
    """
    
    m1 = args['mass1']
    m2 = args['mass2']
    ri = args['ri']
    Qi = args['Qi']
    Dl = args['distance']
    duration = args['duration']
    dt = args['dt']
    vi = np.sqrt(G*(5*10**5)*Ms/(3* 10*pc))

    ### some factors and variables
    Mu, GM, L, phi0, rmin, p, ecc = get_orbit(m1, m2, ri, Qi, vi)
    C1 = (ecc-1)/np.sqrt(ecc**2-1)
    C2 = ecc*np.sqrt(ecc**2-1)/2
    factor = 2*p**2/(ecc**2-1)**(3/2)/L

    # Generate Time array
    times = np.arange(-duration/2, duration/2+dt/2, dt)
    Times = times[:int(len(times)//2)+1]
    
    ## find phit from time using Equaiton function with bisection method
    # set bisection accuracy
        
    Ltf = Times/factor
    phit = [bisect(Equation, 0, np.pi, args=(phi0, ltf, C1, C2, ecc)) for ltf in Ltf]

    phit2 = 2*phi0 - phit[:-1][::-1]
    phit = np.array(list(phit) + list(phit2))
    
    D = 10**Dl * pc
    h11, h22, h12 = H_pols(phit, phi0, Mu, GM, D, L, ecc, rmin)

    return TimeSeries(h11, delta_t=dt) - TimeSeries(h22, delta_t=dt), 2*TimeSeries(h12, delta_t=dt)

def hyperbolic_fd(**args):
    
    hp, hc = hyperbolic_td(**args)
    return hp.to_frequencyseries(), hc.to_frequencyseries()