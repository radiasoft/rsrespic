
import numpy as np
from scipy.constants import c as c_mks

c = c_mks*1.e2

class Field(object):

    def __init__(self, lambda0, numx, numy):
        """Initialize a field object with Fourier basis

        Arguments:
            lambda0 (float) : longest wavelength to be considered in the Fourier basis

        """
        self.lambda0 = lambda0
        self.modes = []
        self.sigma_mn = []
        self.num_x = numx
        self.num_y = numy
        self.frequencies = []

        self.create_modes()

    def create_modes(self, **kwargs):

        """
        Specify mode numbers and instantiate corresponding k values

        Arguments:
            numx (int) : number of modes in x
            numy (int) : number of modes in y
        """

        if kwargs.__contains__('numx'):
            numx = kwargs['numx']
        else:
            numx = self.num_x

        if kwargs.__contains__('numy'):
            numy = kwargs['numy']
        else:
            numy = self.num_y

        # Mode frequencies and wavenumbers
        k0 = 2.*np.pi/ self.lambda0 #fundamental mode

        #ks = np.zeros()
        self.kxs = np.linspace(-numx, numx, numx*2 + 1)*k0
        self.kys = np.linspace(-numy, numy, numy*2 + 1)*k0

        self.sigma_mn = np.zeros((numx,numy)) #array of amplitudes numx rows, numy columns

        self.omegaxs = self.kxs*c
        self.omegays = self.kys*c


    def return_phi(self, x, y):

        """
        Compute phi for a given position, assuming sigma matrix is already known

        Arguments:
            x (float) : x position
            y (float : y position

        """

        k_arg = -1.j *(self.kxs*x + self.kys*y)

        exp = np.exp(k_arg)

        self.phi = np.einsum('ij, ')

        phi = self.phi

        #compute kxx + kyy -> mxn matrix of modes
        kxx = np.einsum('m,n->mn', kxs, np.ones(self.num_y)) * x
        kyy = np.einsum('n,m->mn', kys, np.ones(self.num_x)) * y
        exp_arg = kxx + kyy


        # Field components are exp(i(kxx + kyy))
        phi_modes = np.exp(1.j * exp_arg)

        # now multiply by the sigma matrix, component-wise, and sum
        phi_vals = np.sum(np.einsum('mn,mn->mn', phi, phi_modes))

        #phi_amps = np.einsum('mn,mn->mn', phi, phi_modes)