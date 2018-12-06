

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