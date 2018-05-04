import numpy as np
from constants import cgs_constants
from scipy.special import erf

pi = np.pi

## Convert units to cgs from mks
q = cgs_constants['q']
c = cgs_constants['c'] 

class cartesian_2D:

	def __init__(self, L_x = 1.0, L_y = 1.0, L_x_min = 0.01, L_y_min = 0.01, n_modes_x = 10., n_modes_y = 10.):

		""" Input data from the user for the field solver """
		self.lambda_y_0 = L_x * 2.
		self.lambda_x_0 = L_y * 2.
	
		k_x_min = 2.0 * pi / self.lambda_x_0 
		k_y_min = 2.0 * pi / self.lambda_y_0

		self.k_x_max = 2.0 * pi / L_x_min
		self.k_y_max = 2.0 * pi / L_y_min

		n_modes_x = int(self.k_x_max / k_x_min)
		n_modes_y = int(self.k_y_max / k_y_min)
		
		self.n_modes_y = n_modes_y * 2.
		self.n_modes_x = n_modes_x * 2.

		index_x = np.append(-np.arange(n_modes_x)[::-1] - 1, np.arange(n_modes_x) + 1, axis = 0)
		index_y = np.append(-np.arange(n_modes_y)[::-1] - 1, np.arange(n_modes_y) + 1, axis = 0) 

		""" Compute the k_vectors needed to build the k-matrix"""
		self.k_x_vector = index_x * k_x_min
		self.k_y_vector = index_y * k_y_min

		k_x_matrix_sq,k_y_matrix_sq = np.meshgrid(self.k_x_vector ** 2, self.k_y_vector ** 2)

		k_matrix_sq = k_x_matrix_sq + k_y_matrix_sq

		self.k_sq_inv = 1. / k_matrix_sq 

		self.phi = 0.0
		self.Ax = 0.0 
		self.Ay = 0.0
		self.Az = 0.0

	def register_solver(self,field_solver):

		self.solver = field_solver

		return
		

	def write_output(self, fn):


		return 
		

		
		

