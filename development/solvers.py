import numpy as np
from scipy.constants import c, mu_0, m_e
from scipy.constants import epsilon_0 as e_0
from scipy.constants import elementary_charge as q
from scipy.special import erf
from numpy import exp, sin, einsum



pi = np.pi


class field_solver_2D(object):
	""" Class for computing analytic solutions to fields and potentials 
	for a cylidrical gaussian besm""" 

	def __init__(self, lambda_x_0 = 1.0, lambda_y_0 = 1.0, n_modes_x = 10., n_modes_y = 10.):

		""" Input data from the user for the field solver """
		self.lambda_y_0 = lambda_y_0
		self.lambda_x_0 = lambda_x_0
		self.n_modes_y = n_modes_y * 2
		self.n_modes_x = n_modes_x * 2

		index_x = np.append(-np.arange(n_modes_x)[::-1] - 1, np.arange(n_modes_x) + 1, axis = 0)
		index_y = np.append(-np.arange(n_modes_y)[::-1] - 1, np.arange(n_modes_y) + 1, axis = 0)

		""" Compute the k_vectors needed to build the k-matrix"""
		self.k_x_vector = index_x * 2.0 * pi / self.lambda_x_0
		self.k_y_vector = index_y * 2.0 * pi / self.lambda_y_0

		self.k_x_matrix_sq,self.k_y_matrix_sq = np.meshgrid(self.k_x_vector ** 2, self.k_y_vector ** 2)

		self.k_matrix_sq = self.k_x_matrix_sq + self.k_y_matrix_sq

		self.k_sq_inv = 1. / self.k_matrix_sq # np.linalg.pinv(self.k_matrix_sq)

		#self.k_sq_inv[n_modes_x , n_modes_y ] = 0

	def compute_phi(self, particles):	
		
		kx_phases = einsum('x, p -> xp', self.k_x_vector, particles.x)
		ky_phases = einsum('y, p -> yp', self.k_y_vector, particles.y)

		kx4 = np.einsum('m,n,p -> mnp', self.k_x_vector,np.ones(self.n_modes_y), particles.x)
		ky4 = np.einsum('m,n,p -> mnp', np.ones(self.n_modes_x),self.k_y_vector, particles.y)
		exponential_arg = kx4 + ky4		

		ptcl_exponential = exp(1j * exponential_arg) #* dk  #* self.f_sigma(kx4, ky4, particles)

		phi = einsum('xyp, xy, p -> xy', ptcl_exponential, self.k_sq_inv, particles.charge)

		self.phi = - phi / e_0  #/ (self.lambda_x_0 * self.lambda_y_0)

		return phi


	def compute_phi_mesh(self, **kwargs):
	 	
		if "xmax" in kwargs:
			xmax = kwargs["xmax"]
		else:
		 	xmax = self.lambda_x_0
	
		if "ymax" in kwargs:
			ymax = kwargs["ymax"]
		else:
		 	ymax = self.lambda_y_0
	
		xarray = np.linspace(-xmax, xmax, 100)
		yarray = np.linspace(-ymax, ymax, 100)
		XX, YY = np.meshgrid(xarray, yarray)

		phi = self.phi

		kx4 = np.einsum('m,n,ij -> mnij', self.k_x_vector,np.ones(self.n_modes_y), XX)
		ky4 = np.einsum('m,n,ij -> mnij', np.ones(self.n_modes_x),self.k_y_vector, YY)
		exp_arg = kx4 + ky4

		#Field components are exp(-i(kxx + kyy))
		phi_modes = np.exp(1.j*exp_arg) #+ np.exp(1.j * exp_arg_2)

		#now multiply by the sigma matrix, component-wise - using shared mn indices
		phi_amps = np.einsum('mn,mnij->mnij',phi, phi_modes)

		#now sum over all modes (mn)
		phi_vals = np.einsum('mnij->ij',phi_amps)


		return phi_vals - np.min(phi_vals), XX, YY

	def f_sigma(self, k_x, k_y, particles):

		arg_x = k_x * particles.dx_tent / (2. * pi) 
		arg_y = k_y * particles.dy_tent / (2. * pi)
 
		lambda_twiddle = np.sinc(arg_x) ** 2 * np.sinc(arg_y) ** 2 * particles.dx_tent * particles.dy_tent / (2. * pi)

		return lambda_twiddle






	


