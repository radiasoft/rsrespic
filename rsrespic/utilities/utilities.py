from __future__ import division

## Python imports
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import os
import scipy.integrate as sint
import subprocess

## respic imports 
import constants


q = constants.cgs_constants['q']
c = constants.cgs_constants['c']
m_e = constants.cgs_constants['m_e']
m_p = constants.cgs_constants['m_p']


def compute_system_energy(particles, fields):
    
    h_drift = np.sum(np.sqrt((particles.beta * particles.p_xi)**2 - particles.px ** 2 
                    - particles.py ** 2 - (particles.weight * particles.m_0 * c) ** 2)) * 9.9
    
    fields.solver.compute_psi_particles(fields, particles)
    
    h_sc = -np.sum( particles.weight * particles.charge / (particles.beta * c) * fields.psi_vals)
    
    #es_energy = fields.mode_coefficients * np.conj(fields.mode_coefficients) * fields.k_matrix_sq
    
    return np.abs(h_drift + h_sc) # + np.sum(es_energy))


def plot_beam(beam):

    plt.figure(figsize = (12,12))
    plt.subplot(2,2,1)
    plt.plot(beam.x, beam.px / beam.pz, 'o')
    
    plt.subplot(2,2,2)
    plt.plot(beam.y, beam.py / beam.pz, 'o')
    
    plt.subplot(2,2,3)
    plt.plot(beam.x, beam.y, 'o')
    
    plt.subplot(2,2,4)
    plt.plot(beam.px / beam.pz, beam.py / beam.pz, 'o')
    
    plt.show()
    

def plot_rs_beam(points):
    x = points[:,0]
    xp = points[:,1]
    y = points[:,2]
    yp = points[:,3]


    plt.figure(figsize = (16,16))
    plt.subplot(3,3,1)
    plt.hexbin(x, xp, gridsize = 40)#, 'o')
    plt.title('x-xp')

    plt.subplot(3,3,2)
    plt.hexbin(y, yp, gridsize = 40)#,'o')
    plt.title('y-yp')

    plt.subplot(3,3,3)
    plt.hexbin(x, y, gridsize = 40)#,'o')
    plt.title('x-y')

    plt.subplot(3,3,4)
    plt.hexbin(xp, yp, gridsize = 40)#, 'o')
    plt.title('xp-yp')

    plt.subplot(3,3,5)
    plt.hexbin(x, yp, gridsize = 40)#,'o')
    plt.title('x-yp')

    plt.subplot(3,3,6)
    plt.hexbin(y, xp, gridsize = 40)#, 'o')
    plt.title('y-xp')
    plt.show()

    return



def round_beam_expansion(z, e_x, e_y, x0, y0, q, gamma, m_0):
    
    beta = np.sqrt(1. - 1./ gamma **2)
    
    I = q * beta * c / 100. 
    
    I0 = 4 * np.pi * 8.85e-12 * (c/100.) ** 3 * (m_0 / 1000.) / (1.602e-19) 
    
    K = I * 2. / (I0 * (beta **3) * (gamma **3))

    def func(E, z):
        
        x = E[0]
        xp = E[1]
        y = E[2]
        yp = E[3]
        
        xpp = 2 * K / (x + y) + e_x**2 / x**3
        ypp = 2 * K / (x + y) + e_y**2 / y**3
        
        return np.asarray([xp, xpp, yp, ypp])
    
    init = np.asarray([x0, 0, y0, 0])
    
    output = sint.odeint(func, init, z)
    
    return output



def calc_perveance(I, beam_beta, beam_gamma, cn=0):
    '''Calculate the perveance for a proton beam of a given current and particle energy.
    
    Arguments
        - I - current in A
        - ref - the reference particle for extracting beta and gamma
        
        - (optional) charge neutralization factor - default 0
    '''
    
    I0 = 3.13e7 #characteristic current
    
    beta = beam_beta
    gamma = beam_gamma
    
    return (I/I0)*(2/beta**3)*(1/gamma**3)

def calc_characteristic_current():
    '''Return characteristics current for proton beam'''
    return 4*np.pi*scipy.constants.epsilon_0*scipy.constants.m_p*(scipy.constants.c**3)/scipy.constants.e


#Introduce numerical integrators

#2nd Order RK - Ralston Method
def Ralston(r,z,h,f):
    k1 = h*f(r)
    return 0.25*k1 + 0.75*h*f(r+(2/3)*k1)

#4th Order Runge-Kutta
def RungeKutta4(r,z,h,f):
    k1 = f(r)
    k2 = f(r + (h/2)*k1)
    k3 = f(r + (h/2)*k2)
    k4 = f(r + h*k3)
    return h/6*(k1 + 2*k2 +2*k3 + k4)

#function here, which is a function of r and z
def rprime(K,emit,r0,rp0,rm):
    '''
    
    Returns the slope of the beam envelope (dr/dz) for a given value of emittance,rm, K, and initial conditions.
    
    This equation follows from Reisier.
    
    Arguments:
    
        - r - beam radius (or RMS)
        - K - perveance
        - emit - geometric emittance
        - r0 - initial envelope radius (or RMS)
        - rp0 - initial slope of envelope (or RMS)
        
    '''
    
    first = rp0**2 #first term
    second = (emit**2)*((1./r0**2)-(1./rm**2)) #second term
    third = 2*K* np.log(rm/r0) / 4
    
    return np.sqrt(first + second + third)


def calculate_expansion(current, beam_beta, beam_gamma, r0, rprime0, 
    emit = 1.0e-6, N = 1000, zf = 1.0):

    '''Evaluate the expansion of a KV beam envelope in a drift along z-axis, begining at z = 0.
    
    Arguments:
        - current - beam current in A
        - reference_particle - synergia object for bunch/lattice reference particle
        - r0 - initial envelope value (provide RMS for RMS expansion, a for envelope expansion, etc.)
        - rp0 - initial slope of envelope (must be non-zero, but calculation is not sensitive to small values)
        
        - (optional) emit - geometric emittance of beam - default 2.05721258396*1.e-6 (for 0.3 mm-mrad KV beam)
        - (optional) N - number of steps for integration - default 1000
        - (optional) zf - final z value (e.g. length of expansion) - default 50.0
        
    '''
    
    z0 = 0.0 #start
    ss = (zf-z0)/N #step size

    zpoints = np.linspace(0.0, zf, num=N) #define z values
    rpoints = [] #empty array for r values
    
    #calculate perveance
    Kp = calc_perveance(current, beam_beta, beam_gamma)
    
    #x is r
    #z is t (what we step up)
    #f is our function describing the relationship between r and z
    f = lambda r: rprime(Kp,emit,r0,rprime0,r)

    r,z,dz = r0,z0,ss
    points = []
    while z < zf:
        points.append((z,r))
        z, r = z+dz, r + Ralston(r,z,dz,f) #incremement
        
    return points



def read_sdds_columns(file_name, column_names, data_class = None):

    if data_class == None:
        class sdds_data:
            def __init__(self,file_name,column_names):
                self.file_name = file_name
                self.column_names = column_names

        data_class = sdds_data(file_name,column_names)


    for i in range(0,len(column_names)):
        column = column_names[i]
        stream = subprocess.Popen("""sdds2stream """+file_name+""" -col="""+column+""" -delimiter=" "
        """ , shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)   
        

        column_data = np.asarray(stream.stdout.read().split())
        error = stream.stderr.read()

        setattr(data_class, column, column_data.astype(float))
        setattr(data_class, column + '.error', error)

    return data_class



def plot_twiss_sdds(beam, flag = None):
    
    if flag == None:
        plt.figure(figsize = (12,12))
        plt.subplot(2,1,1)
        plt.plot(beam.s, beam.alphax)
        plt.plot(beam.s, beam.alphay)
        plt.legend(['alphax', 'alphay'])
        plt.xlabel('s')

        plt.subplot(2,1,2)
        plt.plot(beam.s, beam.betax)
        plt.plot(beam.s, beam.betay)
        plt.legend(['betax', 'betay'])
        plt.xlabel('s')

        plt.show()

    if flag == 'Beam':
        plt.figure(figsize = (12,12))
        plt.subplot(2,1,1)
        plt.plot(beam.s, beam.alphaxBeam)
        plt.plot(beam.s, beam.alphayBeam)
        plt.legend(['alphax', 'alphay'])
        plt.xlabel('s')

        plt.subplot(2,1,2)
        plt.plot(beam.s, beam.betaxBeam)
        plt.plot(beam.s, beam.betayBeam)
        plt.legend(['betax', 'betay'])
        plt.xlabel('s')

        plt.show()


def plot_elegant_beam(beam):

    plt.figure(figsize = (16,16))
    plt.subplot(3,3,1)
    plt.hexbin(beam.x, beam.xp, gridsize = 40)#, 'o')
    plt.title('x-xp')

    plt.subplot(3,3,2)
    plt.hexbin(beam.y, beam.yp, gridsize = 40)#,'o')
    plt.title('y-yp')

    plt.subplot(3,3,3)
    plt.hexbin(beam.x, beam.y, gridsize = 40)#,'o')
    plt.title('x-y')

    plt.subplot(3,3,4)
    plt.hexbin(beam.xp, beam.yp, gridsize = 40)#, 'o')
    plt.title('xp-yp')

    plt.subplot(3,3,5)
    plt.hexbin(beam.x, beam.yp, gridsize = 40)#,'o')
    plt.title('x-yp')

    plt.subplot(3,3,6)
    plt.hexbin(beam.y, beam.xp, gridsize = 40)#, 'o')
    plt.title('y-xp')
    plt.show()

    return
