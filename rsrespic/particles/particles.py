# Particle class for rsrespic
#
# Class is initialized with a bunch of particles to provide position and momentum values
# Weightings are optional as the initial implementation is gridless
#
# Keep in mind that s is the independent variable, and thus x,y, and tau = c*t are the coordinate
# descriptors of the particles, along with px, py, and ptau.
#
# Coordinates:
#   - Particles have position x, y, tau
#
# Momenta
#   - Particles have momenta p_x p_y, p_tau
#   - We assume a fixed beta for the bunch, computed by averaging over all tau velocities
#
# Units
#   -Assume CGS units for now
#

import numpy as np
from scipy.constants import m_e as me_mks
from scipy.constants import c as c_mks

# Set the default mass and charge for an electron
m = me_mks*1.e3 #cgs
q = 4.80320451e-10 #esu 1.*e
c = c_mks*1.e2 #cgs


class Particles(object):

    def __init__(self, n_particles, charge, mass, weight, **kwargs):

        """
        Class that stores the particle data and can compute the particle ensemble quantities

        Arguments
        ----------
            n_total (int) : total number of macroparticles across all processors

            n_particles (int) : number of macroparticles on local processor

            charge float (float) : charge of the physical species in units e

            mass: float (grams) : mass of the physical species

            weight float : number of particles per macroparticle

            shape: string (optional): specifying the particle shape function: defaults to tent function for now

            species_name string (optional) name of the particle species
        """


        self.n_total = n_particles
        self.np = n_particles

        if kwargs.__contains__('n_total'):
            self.n_total = kwargs['n_total']

        if kwargs.__contains__('species_name'):
            self.species_name = kwargs['species_name']

        #dynamical quantities
        self.x = np.zeros(n_particles)
        self.y = np.zeros(n_particles)
        self.z = np.zeros(n_particles)
        self.px = np.zeros(n_particles)
        self.py = np.zeros(n_particles)
        self.pz = np.zeros(n_particles)

        self.gamma = np.zeros(n_particles)
        self.gamma_mc = np.zeros(n_particles)

        self.weight = weight * np.ones(n_particles)

        self.charge = charge*q
        self.mass = mass

        self.qOc = (self.charge / c) * np.ones(n_particles)
        self.mc = (self.mass * c) * np.ones(n_particles)


        if kwargs.__contains__('shape'):
            self.shape = kwargs['shape']
        else:
            self.shape = 'delta' #'tent'


    def add_particles(self, positions, momenta, weights=None, IDs=None):
        '''
        Initialize bunch of particles. Overwrite position and momentum arrays

        Arguments:
            positions (ndarray): array of positions - [x, y, z]

            momenta (ndarray): array of momenta - [px, py, pz]

            weights (Optional[ndarray]): array of weights- [wx,wy,wz]

            IDs (Optional[ndarray]): array of particle IDs - length # of particles

        '''

        if not type(positions) == 'ndarray':
            positions = np.asarray(positions)
        if not type(momenta) == 'ndarray':
            momenta = np.asarray(momenta)
        if not positions.shape[0] == momenta.shape[0]:
            print "Position and momentum arrays have unequal lengths"
            raise

        # Position quantities
        self.x = positions[:, 0]
        self.y = positions[:, 1]
        self.z = positions[:, 2]

        self.num_particles = len(self.x)

        # initialize particle IDs
        if not IDs is None:
            if len(IDs) == self.num_particles:
                self.IDs = IDs
            else:
                print "Number of particle IDs differs from number of particles"
                raise
        else:
            self.IDs = np.arange(self.num_particles) + 1

        # initialize weights
        if weights is None:
            self.weights = np.ones(self.num_particles)
        elif not type(weights) == 'ndarray':
            weights = np.asarray(weights)
            if len(weights) == self.num_particles:
                self.weights = weights
            else:
                print "Number of particle weights differs from number of particles"
                raise

        # Charge and mass quantities - weighted
        #self.mass = self.weights * self.mass
        self.qOc = self.weights * self.q0c
        self.mc = self.weights * self.mc

        # Momentum quantities - weighted
        self.px = self.weights * momenta[:, 0]
        self.py = self.weights * momenta[:, 1]
        self.pz = self.weights * momenta[:, 2]

    def construct_gaussian(self, sigma_x=1.0, sigma_y=1.0):

        sigma_x = sigma_x
        sigma_y = sigma_y
        sigma_xp = 0.01
        sigma_yp = 0.01

        sigma_xxp = 0.0005
        sigma_xy = 0.0000
        sigma_xyp = 0.0000
        sigma_xpy = 0.0000
        sigma_xpyp = 0.0000
        sigma_yyp = 0.0000

        mean = [0, 0, 0, 0]

        cov = [[sigma_x ** 2, sigma_xxp ** 2, sigma_xy ** 2, sigma_xyp ** 2],
               [sigma_xxp ** 2, sigma_xp ** 2, sigma_xpy ** 2, sigma_xpyp ** 2],
               [sigma_xy ** 2, sigma_xpy ** 2, sigma_y ** 2, sigma_yyp ** 2],
               [sigma_xyp ** 2, sigma_xpyp ** 2, sigma_yyp ** 2, sigma_yp ** 2]]

        x, px, y, py = np.random.multivariate_normal(mean, cov, self.np).T

        self.x = x
        self.y = y
        #self.charge = self.w * np.ones(len(x))