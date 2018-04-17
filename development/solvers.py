import numpy as np
from scipy.constants import c, mu_0, m_e
from scipy.constants import epsilon_0 as e_0
from scipy.constants import elementary_charge as q
from scipy.special import erf
from numpy import exp, sin, einsum

import colormaps as cmaps
import matplotlib.pyplot as plt

pi = np.pi


class field_solver_2D(object):
	""" Class for computing analytic solutions to fields and potentials 
	for a cylidrical gaussian besm""" 

	def __init__(self, L_x = 1.0, L_y = 1.0, L_x_min = 0.01, L_y_min = 0.01, n_modes_x = 10., n_modes_y = 10.):

		""" Input data from the user for the field solver """
		self.lambda_y_0 = L_x * 2.
		self.lambda_x_0 = L_y * 2.
	

		k_x_min = 2.0 * pi / self.lambda_x_0 
		k_y_min = 2.0 * pi / self.lambda_y_0

		k_x_max = 2.0 * pi / L_x_min
		k_y_max = 2.0 * pi / L_y_min

		n_modes_x = int(k_x_max / k_x_min)
		n_modes_y = int(k_y_max / k_y_min)

		print n_modes_y, n_modes_x
		
		self.n_modes_y = n_modes_y * 2.
		self.n_modes_x = n_modes_x * 2.

		index_x = np.append(-np.arange(n_modes_x)[::-1] - 1, np.arange(n_modes_x) + 1, axis = 0)
		index_y = np.append(-np.arange(n_modes_y)[::-1] - 1, np.arange(n_modes_y) + 1, axis = 0) 

		""" Compute the k_vectors needed to build the k-matrix"""
		self.k_x_vector = index_x * k_x_min
		self.k_y_vector = index_y * k_y_min

		k_x_matrix_sq,k_y_matrix_sq = np.meshgrid(self.k_x_vector ** 2, self.k_y_vector ** 2)

		k_matrix_sq = k_x_matrix_sq + k_y_matrix_sq

		self.k_sq_inv = 1. / k_matrix_sq # np.linalg.pinv(self.k_matrix_sq)

		#self.k_sq_inv[n_modes_x , n_modes_y ] = 0

	def compute_phi(self, particles):	
		
		kx_phases = einsum('x, p -> xp', self.k_x_vector, particles.x)
		ky_phases = einsum('y, p -> yp', self.k_y_vector, particles.y)

		kx1 = np.einsum('m, p -> mp', self.k_x_vector, particles.x) 
		ky1 = np.einsum('n, p -> np', self.k_y_vector, particles.y) 

		exp_x = np.exp(1j * kx1)
		exp_y = np.exp(1j * ky1)
	
		ptcl_exponential = np.einsum('mp, np -> mn', exp_x, exp_y) * particles.w #* self.f_sigma(kx4, ky4, particles)

		phi = einsum('xy, xy -> xy', ptcl_exponential, self.k_sq_inv) * 8. * pi

		self.phi = - phi / e_0  

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

		if "n_grid" in kwargs:
			n_grid = kwargs["n_grid"]
		else:
			n_grid = 10


		xarray = np.linspace(-xmax, xmax, n_grid)
		yarray = np.linspace(-ymax, ymax, n_grid)
		XX, YY = np.meshgrid(xarray, yarray)

		phi = self.phi / (self.lambda_y_0 * self.lambda_x_0)

		kx4 = np.einsum('m,i -> mi', self.k_x_vector, xarray)
		ky4 = np.einsum('n,j -> nj', self.k_y_vector, yarray)

		exp_x = np.exp(-1j * kx4)
		exp_y = np.exp(-1j * ky4)

		#Field components are exp(-i(kxx + kyy))
		phi_modes = np.einsum('mi, nj -> mnij', exp_x, exp_y)

		#now multiply by the sigma matrix, component-wise - using shared mn indices
		phi_vals = np.einsum('mn,mnij->ij',phi, phi_modes)


		return phi_vals - np.min(phi_vals), XX, YY

	def f_sigma(self, k_x, k_y, particles):

		arg_x = k_x * particles.dx_tent / (2. * pi) 
		arg_y = k_y * particles.dy_tent / (2. * pi)
 
		lambda_twiddle = np.sinc(arg_x) ** 2 * np.sinc(arg_y) ** 2 * particles.dx_tent * particles.dy_tent / (2. * pi)

		return lambda_twiddle






	


