
import numpy as np


class uniform_focusing_channel:

    def __init__(self, s0 = 0, k = 0, L = 0):

        self.s0 = s0
        self.k = k
        self.L = L

    def thin_quad(self, particles, ds):
        kappa = self.k

        ## Thin quad has no drift 
        dx_ds = 0.
        dy_ds = 0.
        dz_ds = 0.

        ## assume a positron convention here kappa is the same as elegant K1
        dpx = -kappa * particles.x * ds * (particles.pz / 10000.)
        dpy = -kappa * particles.y * ds * (particles.pz / 10000.)
        dpz = 0.

        particles.px += dpx
        particles.py += dpy
        particles.pt += 0. 

    def one_turn_map(self, maps, fields, particles, diagnostics, s, dumper = None):
       
        maps.drift(particles, ds = self.L/2.)
        self.thin_quad(particles, ds = self.L/2.)
        maps.space_charge_kick_2D_sine(fields, particles, ds = self.L)
        self.thin_quad(particles, ds = self.L/2.)
        maps.drift(particles, ds = self.L/2.)
        s += self.L
        diagnostics.update(s, particles)

        return s

    def one_turn_map_cartesian(self, maps, fields, particles, diagnostics, s, dumper = None):
       
        maps.drift(particles, ds = self.L/2.)
        self.thin_quad(particles, ds = self.L/2.)
        maps.space_charge_kick_2D(fields, particles, ds = self.L)
        self.thin_quad(particles, ds = self.L/2.)
        maps.drift(particles, ds = self.L/2.)
        s += self.L
        diagnostics.update(s, particles)

        return s

    
    

class octo_respic: 
    def __init__(self, s0 = 0):
        ## Drit defined in meters and converted to CM for respic
        self.s0 = s0


    def drift(self, maps, fields, particles, diagnostics, s, L = 0):
    
        maps.drift(particles, ds = L/2.)

        #sigma_x = np.std(particles.x)
        #sigma_y = np.std(particles.y)
        #L_x = 30. * sigma_x 
        #L_y = 30. * sigma_y
        #fields.reset_modes(L_x = L_x, L_y = L_y,
        #L_x_min = L_x/30., L_y_min = L_y/30.)
        
        maps.space_charge_kick_2D(fields, particles, ds = L)
        maps.drift(particles, ds = L/2.)
        s += L
        diagnostics.update(s, particles)
    
        return s

    def quad(self, maps, fields, particles, diagnostics, s, L = 0., k = 0):

        s = self.drift(maps, fields, particles, diagnostics, s, L = L/4.)
        maps.thin_quad(fields, particles, ds = L/2., kappa = k)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L/4.)

        s = self.drift(maps, fields, particles, diagnostics, s, L = L/4.)
        maps.thin_quad(fields, particles, ds = L/2., kappa = k)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L/4.)

        #s += L * 100.
        #diagnostics.update(s, particles) 
        
        return s

    def drift_1(self, maps, fields, particles, diagnostics, s, dumper = None):
        L_d1 = 0.5 # units of meters

        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
    
        return s


    def one_turn_map(self, maps, fields, particles, diagnostics, s, dumper = None):
        L_quad = 0.1

        s = self.quad(maps, fields, particles, diagnostics, s, L = L_quad / 2. * 100., k = 0.5)

        s = self.drift_1(maps, fields, particles, diagnostics, s, dumper = dumper)

        s = self.quad(maps, fields, particles, diagnostics, s, L = L_quad * 100., k = -0.5)
        
        s = self.drift_1(maps, fields, particles, diagnostics, s, dumper = dumper)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = L_quad / 2.* 100., k = 0.5)

        return s


class fodo_respic_syn: 
    def __init__(self, s0 = 0):
        ## Drit defined in meters and converted to CM for respic
        self.s0 = s0


    def drift(self, maps, fields, particles, diagnostics, s, L = 0):
    
        maps.drift(particles, ds = L/2.)
        maps.space_charge_kick_2D(fields, particles, ds = L)
        maps.drift(particles, ds = L/2.)
        s += L
        diagnostics.update(s, particles)
    
        return s

    def quad(self, maps, fields, particles, diagnostics, s, L = 0., k = 0):

        maps.thin_quad(fields, particles, ds = L/4., kappa = k)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L/4.)
        maps.thin_quad(fields, particles, ds = L/4., kappa = k)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L/4.)
        maps.thin_quad(fields, particles, ds = L/4., kappa = k)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L/4.)
        maps.thin_quad(fields, particles, ds = L/4., kappa = k)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L/4.)

        #s += L * 100.
        #diagnostics.update(s, particles) 
        
        return s

    def drift_1(self, maps, fields, particles, diagnostics, s):
        L_d1 = 0.5 # units of meters
        L_d2 = 1.0 # units of meters

        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d2 / 2. * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d2 / 2. * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = L_d1 * 100.)
        
        return s


    def one_turn_map(self, maps, fields, particles, diagnostics, s):
        L_quad = 0.1

        s = self.quad(maps, fields, particles, diagnostics, s, L = L_quad / 2. * 100., k = 3.5)

        s = self.drift_1(maps, fields, particles, diagnostics, s)

        s = self.quad(maps, fields, particles, diagnostics, s, L = L_quad * 100., k = -4.0)
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = L_quad / 2.* 100., k = 3.5)

        return s



class fodo_cell:

    def __init__(self, L = 0.1, s0 = 0, k_array = [0, 0, 0, 0]):
        ## Drit defined in meters and converted to CM for respic
        self.s0 = s0
        self.L = L 
        self.k_array = k_array


    def drift(self, maps, fields, particles, diagnostics, s, L = 0):
    
        maps.drift(particles, ds = L)
        s += L
        diagnostics.update(s, particles)
    
        return s

    def quad(self, maps, fields, particles, diagnostics, s, L = 1.0e-6, k = 0):

        maps.thin_quad(fields, particles, ds = L * 100., kappa = k)
        s += L * 100.
        diagnostics.update(s, particles) 
        
        return s

    def drift_1(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        
        return s

    def drift_2(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)

        return s


    def one_turn_map(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)

        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[0])
        
        s = self.drift_2(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[1])

        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[2])

        s = self.drift_2(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[3])

        s = self.drift_1(maps, fields, particles, diagnostics, s)

        return s


class fodo_cell2:

    def __init__(self, L = 0.1, s0 = 0, k_array = [0, 0, 0, 0]):
        ## Drit defined in meters and converted to CM for respic
        self.s0 = s0
        self.L = L 
        self.k_array = k_array


    def drift(self, maps, fields, particles, diagnostics, s, L = 0):
    
        maps.drift(particles, ds = L)
        s += L
        diagnostics.update(s, particles)
    
        return s

    def quad(self, maps, fields, particles, diagnostics, s, L = 1.0e-6, k = 0):

        maps.thin_quad(fields, particles, ds = L * 100., kappa = k)
        s += L * 100.
        diagnostics.update(s, particles) 
        
        return s

    def drift_1(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        
        return s


    def one_turn_map(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)

        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[0])
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[1])

        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[2])

        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[3])

        s = self.drift_1(maps, fields, particles, diagnostics, s)

        return s


class fodo_cell_2:

    def __init__(self, L = 0.1, s0 = 0, k_array = [0, 0, 0, 0]):
        ## Drit defined in meters and converted to CM for respic
        self.s0 = s0
        self.L = L 
        self.k_array = k_array


    def drift(self, maps, fields, particles, diagnostics, s, L = 0):
    
        maps.drift(particles, ds = L)
        s += L
        diagnostics.update(s, particles)
    
        return s

    def quad(self, maps, fields, particles, diagnostics, s, L = 1.0e-6, k = 0):

        maps.thin_quad(fields, particles, ds = L * 100., kappa = k)
        s += L * 100.
        diagnostics.update(s, particles) 
        
        return s

    def drift_1(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        
        return s

    def one_turn_map(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)

        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[0])
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[1])

        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[2])

        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[3])

        s = self.drift_1(maps, fields, particles, diagnostics, s)

        return s



class fodo_cell_sc:

    def __init__(self, L = 0.1, s0 = 0, k_array = [0, 0, 0, 0]):
        ## Drit defined in meters and converted to CM for respic
        self.s0 = s0
        self.L = L 
        self.k_array = k_array


    def drift(self, maps, fields, particles, diagnostics, s, L = 0):
        
        maps.drift(particles, ds = L/2.)
        maps.space_charge_kick_2D(fields, particles, ds = L)
        maps.drift(particles, ds = L/2.)
        s += L
        diagnostics.update(s, particles)
    
        return s

    def quad(self, maps, fields, particles, diagnostics, s, L = 1.0e-6, k = 0):

        maps.thin_quad(fields, particles, ds = L * 100., kappa = k)
        s += L * 100.
        diagnostics.update(s, particles) 
        
        return s

    def drift_1(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        
        return s

    def drift_2(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)

        return s


    def one_turn_map(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)

        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[0])
        
        s = self.drift_2(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[1])

        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[2])

        s = self.drift_2(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[3])

        s = self.drift_1(maps, fields, particles, diagnostics, s)

        return s


class fodo_cell2_sc:

    def __init__(self, L = 0.1, s0 = 0, k_array = [0, 0, 0, 0]):
        ## Drit defined in meters and converted to CM for respic
        self.s0 = s0
        self.L = L 
        self.k_array = k_array


    def drift(self, maps, fields, particles, diagnostics, s, L = 0):
        
        maps.drift(particles, ds = L/2.)
        maps.space_charge_kick_2D(fields, particles, ds = L)
        maps.drift(particles, ds = L/2.)
        s += L
        diagnostics.update(s, particles)
    
        return s

    def quad(self, maps, fields, particles, diagnostics, s, L = 1.0e-6, k = 0):

        maps.thin_quad(fields, particles, ds = L * 100., kappa = k)
        s += L * 100.
        diagnostics.update(s, particles) 
        
        return s

    def drift_1(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        s = self.drift(maps, fields, particles, diagnostics, s, L = self.L * 100.)
        
        return s


    def one_turn_map(self, maps, fields, particles, diagnostics, s):
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)

        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[0])
        
        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[1])

        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[2])

        s = self.drift_1(maps, fields, particles, diagnostics, s)
        
        s = self.quad(maps, fields, particles, diagnostics, s, L = 1.0e-6, k = self.k_array[3])

        s = self.drift_1(maps, fields, particles, diagnostics, s)

        return s