
import numpy as np
import matplotlib
import scipy.integrate as sint
import matplotlib.pyplot as plt

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