
import numpy as np
from scipy.constants import c as c_mks

c = c_mks*1.e2
pi = np.pi

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


class sin_transform_2D:

    def __init__(self, L_x = 1.0, L_y = 1.0, n_modes_x = 10., n_modes_y = 10.):

        """ Input data from the user for the field solver """

        self.L_x = L_y * 2.
        self.L_y = L_x * 2.
        self.n_modes_x = n_modes_x
        self.n_modes_y = n_modes_y

        self.phi = 0.0
        self.Ax = 0.0 
        self.Ay = 0.0
        self.Az = 0.0

    def register_solver(self,field_solver):

        self.solver = field_solver

        return


class cartesian_2D:

    def __init__(self, L_x = 1.0, L_y = 1.0, L_x_min = 0.01, L_y_min = 0.01, n_modes_x = 10., n_modes_y = 10.):

        """ Input data from the user for the field solver """
        self.lambda_y_0 = L_y * 2.
        self.lambda_x_0 = L_x * 2.
    
        k_x_min = 2.0 * pi / self.lambda_x_0 
        k_y_min = 2.0 * pi / self.lambda_y_0

        self.k_x_max = 2.0 * pi / L_x_min
        self.k_y_max = 2.0 * pi / L_y_min

        n_modes_x = np.rint(self.k_x_max / k_x_min)
        n_modes_y = np.rint(self.k_y_max / k_y_min)
        
        self.n_modes_y = n_modes_y * 2.
        self.n_modes_x = n_modes_x * 2.

        index_x = np.append(-np.arange(n_modes_x)[::-1] - 1, np.arange(n_modes_x) + 1, axis = 0)
        index_y = np.append(-np.arange(n_modes_y)[::-1] - 1, np.arange(n_modes_y) + 1, axis = 0) 

        """ Compute the k_vectors needed to build the k-matrix"""
        self.k_x_vector = index_x * k_x_min
        self.k_y_vector = index_y * k_y_min

        k_x_matrix_sq,k_y_matrix_sq = np.meshgrid(self.k_x_vector ** 2, self.k_y_vector ** 2)

        k_matrix_sq = k_x_matrix_sq + k_y_matrix_sq

        self.k_matrix_sq = k_matrix_sq

        self.k_sq_inv = 1. / k_matrix_sq 

        self.phi = 0.0
        self.Ax = 0.0 
        self.Ay = 0.0
        self.Az = 0.0

    def register_solver(self,field_solver):

        self.solver = field_solver

        return
        

    def reset_modes(self, L_x = 1.0, L_y = 1.0, L_x_min = 0.01, L_y_min = 0.01):

        """ Input data from the user for the field solver """
        self.lambda_y_0 = L_y * 2.
        self.lambda_x_0 = L_x * 2.
    
        k_x_min = 2.0 * pi / self.lambda_x_0 
        k_y_min = 2.0 * pi / self.lambda_y_0

        self.k_x_max = 2.0 * pi / L_x_min
        self.k_y_max = 2.0 * pi / L_y_min

        n_modes_x = np.rint(self.k_x_max / k_x_min)
        n_modes_y = np.rint(self.k_y_max / k_y_min)
        
        self.n_modes_y = n_modes_y * 2.
        self.n_modes_x = n_modes_x * 2.

        index_x = np.append(-np.arange(n_modes_x)[::-1] - 1, np.arange(n_modes_x) + 1, axis = 0)
        index_y = np.append(-np.arange(n_modes_y)[::-1] - 1, np.arange(n_modes_y) + 1, axis = 0) 

        """ Compute the k_vectors needed to build the k-matrix"""
        self.k_x_vector = index_x * k_x_min
        self.k_y_vector = index_y * k_y_min

        k_x_matrix_sq,k_y_matrix_sq = np.meshgrid(self.k_x_vector ** 2, self.k_y_vector ** 2)

        k_matrix_sq = k_x_matrix_sq + k_y_matrix_sq

        self.k_matrix_sq = k_matrix_sq

        self.k_sq_inv = 1. / k_matrix_sq 

        self.phi = 0.0
        self.Ax = 0.0 
        self.Ay = 0.0
        self.Az = 0.0

    def write_output(self, fn):


        return
        