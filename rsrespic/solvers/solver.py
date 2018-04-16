import numpy as np


class ES2D(object):
    
    def __init__(self,fields):
        '''
        Class for solving for the electrostatic potential and its gradient.
        
        Arguments:
            fields instance
            ptcls instance
        
        '''
        
        self.kx = fields.kxs
        self.ky = fields.kys
        self.l0 = fields.lambda0

        # needs to be a good way to calculate this cleanly, otherwise a for loop
        #self.k_sqrd_inv = np.zeros(np.shape(self.kx)[0], np.shape(ky)[0])


        self.k_x_sq, self.k_y_sq = np.meshgrid(self.kx ** 2, self.ky ** 2)

        self.k_matrix_sq = self.k_x_sq + self.k_y_sq

        self.k_sq_inv = 1. / self.k_matrix_sq
        #set the km,n,=0 term to 0 rather than 1/0



        self.k_sq_inv[fields.num_x, fields.num_y] = 0


    def compute_phi(self, ptcls):
        '''Compute the electrostatic potential in terms of Fourier components phi_{ij}
        for a given particle distribution
        '''

        # compute relevant parameters for the algorithm
        #beta0 = ptcls.beta
        #gamma0sqrd = 1/(1. - beta0*beta0)


        #Need to compute the kx*theta terms
        kx_phase = np.einsum('x, p -> xp', self.kx, ptcls.x)
        kx_phase = np.einsum('y, p -> yp', self.ky, ptcls.y)
        
        #kx4/ky4 computes kxx for all particles and all mode combinations
        kx4 = np.einsum('m,n,p -> mnp', self.kx, np.ones(len(self.ky)), ptcls.x)
        ky4 = np.einsum('m,n,p -> mnp', np.ones(len(self.kx)), self.ky, ptcls.y)

        exponential_arg = kx4 + ky4

        #exponential is the sum of both
        ptcl_exponential = np.exp(1j * exponential_arg) #* self.f_sigma(kx4, ky4, ptcls)


        #phi is sum over modes of ksquard*theconvolution
        phi = np.einsum('xyp, xy, p -> xy', ptcl_exponential, self.k_sq_inv, ptcls.qOc)

        self.phi = phi

        return phi



    def compute_phi_mesh(self, **kwargs):
        '''
        Compute potential and apply it to gridded data.
        '''

        if kwargs.__contains__('xmax'):
            xmax = kwargs['xmax']
        else:
            xmax = self.l0

        if kwargs.__contains__('ymax'):
            ymax = kwargs['ymax']
        else:
            ymax = self.l0

        #define number of nodes for applying phi to mesh
        xnodes = 100
        ynodes = 100

        xarray = np.linspace(-xmax, xmax, xnodes)
        yarray = np.linspace(-ymax, ymax, ynodes)
        XX, YY = np.meshgrid(xarray, yarray)

        phi = self.phi

        #permute kxx and kyy over all m,n,i,j indices
        kx4 = np.einsum('m,n,ij -> mnij', self.kx, np.ones(len(self.ky)), XX)
        ky4 = np.einsum('m,n,ij -> mnij', np.ones(len(self.kx)), self.ky, YY)
        exp_arg = kx4 + ky4

        # Field components are exp(-i(kxx + kyy))
        phi_modes = np.exp(1.j * exp_arg)  # + np.exp(1.j * exp_arg_2)

        # now multiply by the sigma matrix, component-wise - using shared mn indices
        phi_amps = np.einsum('mn,mnij->mnij', phi, phi_modes)

        # now sum over all modes (mn)
        phi_vals = np.einsum('mnij->ij', phi_amps)

        return phi_vals, XX, YY

    def f_sigma(self, k_x, k_y, particles):
        '''f sigma computes the spectral composition of the particle shape function (for a tent function)'''

        arg_x = k_x * particles.dx_tent / (2. * pi)
        arg_y = k_y * particles.dy_tent / (2. * pi)

        lambda_twiddle = np.sinc(arg_x) ** 2 * np.sinc(arg_y) ** 2 * particles.dx_tent * particles.dy_tent / (2. * pi)

        return lambda_twiddle
        
        