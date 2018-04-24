import numpy as np
from scipy.constants import c, mu_0, m_e
from scipy.constants import epsilon_0 as e_0
from scipy.constants import elementary_charge as q
from scipy.special import erf
from numpy import exp, sin, einsum

pi = np.pi


class field_solver_2D(object):

	def __init__(self):

		self.name = '2-d electrostatic field solver'

	def compute_phi(self, fields, particles):
		
		#kx_phases = einsum('x, p -> xp', self.k_x_vector, particles.x)
		#ky_phases = einsum('y, p -> yp', self.k_y_vector, particles.y)

		kx1 = np.einsum('m, p -> mp', fields.k_x_vector, particles.y) 
		ky1 = np.einsum('n, p -> np', fields.k_y_vector, particles.x) 

		#kx_mat = np.einsum('m, p -> mp', fields.k_x_vector, np.ones(len(particles.x))) 
		#ky_mat = np.einsum('n, p -> np', fields.k_y_vector, np.ones(len(particles.y)))

		trash, kx_mat = np.meshgrid(particles.x, fields.k_x_vector)
		trash, ky_mat = np.meshgrid(particles.y, fields.k_y_vector)

		exp_x = np.exp(1j * kx1) * particles.lambda_twiddle(kx_mat, particles.x_extent)
		exp_y = np.exp(1j * ky1) * particles.lambda_twiddle(ky_mat, particles.y_extent)
	
		ptcl_exponential = np.einsum('mp, np -> mn', exp_x, exp_y) 

		phi = einsum('xy, xy -> xy', ptcl_exponential, fields.k_sq_inv)

		fields.phi = - phi * particles.charge * particles.weight  * 4. * pi / ( c * (particles.gamma ** 2. - 1) )

		return 


	def compute_phi_mesh(self, fields,  **kwargs):
	 	
		if "xmax" in kwargs:
			xmax = kwargs["xmax"]
		else:
		 	xmax = fields.lambda_x_0
	
		if "ymax" in kwargs:
			ymax = kwargs["ymax"]
		else:
		 	ymax = fields.lambda_y_0

		if "n_grid" in kwargs:
			n_grid = kwargs["n_grid"]
		else:
			n_grid = 10


		xarray = np.linspace(-xmax, xmax, n_grid)
		yarray = np.linspace(-ymax, ymax, n_grid)
		XX, YY = np.meshgrid(xarray, yarray)

		phi = fields.phi / (fields.lambda_y_0 * fields.lambda_x_0)

		kx4 = np.einsum('m,i -> mi', fields.k_x_vector, xarray)
		ky4 = np.einsum('n,j -> nj', fields.k_y_vector, yarray)

		exp_x = np.exp(-1j * kx4)
		exp_y = np.exp(-1j * ky4)

		#Field components are exp(-i(kxx + kyy))
		phi_modes = np.einsum('mi, nj -> mnij', exp_x, exp_y)

		#now multiply by the sigma matrix, component-wise - using shared mn indices
		phi_vals = np.einsum('mn,mnij->ij',phi, phi_modes)

		fields.phi_grid = phi_vals - np.min(phi_vals)
		fields.x_grid = XX
		fields.y_grid = YY

		return


	def compute_kick(self, fields, particles):

		phi = fields.phi  #/ (fields.lambda_y_0 * fields.lambda_x_0)

		kx4 = np.einsum('m,i -> mi', fields.k_x_vector, particles.x)
		ky4 = np.einsum('n,i -> ni', fields.k_y_vector, particles.y)

		exp_x = np.exp(-1j * kx4)		
		exp_y = np.exp(-1j * ky4)

		kick_modes = np.einsum('mp, np -> mnp', exp_x, exp_y)

		kick_x = np.einsum('mn, m, mnp -> p', phi, 1j * fields.k_x_vector, kick_modes)
		kick_y = np.einsum('mn, n, mnp -> p', phi, 1j * fields.k_y_vector, kick_modes)

		fields.kick_x = np.real(kick_x)
		fields.kick_y = np.real(kick_y)

		return


class kinetics_solver_SC2D:

	def __init__(self, ds = 1.0e-2):

		self.ds = ds



	def compute_momentum_arg(self,particles):

		## Calculate the square root term in the hamaltonian
		argumnent = np.sqrt(particles.beta**2 * particles.p_xi **2 - particles.px **2 - particles.py**2 - particles.m_0**2 * c**2)
		
		return argumnent


	def step(self, particles, fields):

		## The first thing is to trasform the longitudional momentum into canaonical coordiantes
		particles.compute_p_xi

		## This computes the square root term of the hamaltonian that is used in all drifts
		argumnent = 1. / self.compute_momentum_arg(particles)
		relativistic_factor = 1. / particles.beta ** 2 - 1. 


		## First half drift
		particles.x += - particles.px * relativistic_factor * argumnent * self.ds / 2.
		particles.y += - particles.py * relativistic_factor * argumnent * self.ds / 2.

		## compute the full kick
		kick_x =  - fields.kick_x  * particles.charge * particles.weight / c
		kick_y =  - fields.kick_y * particles.mass * particles.weight / c

		## apply the kick 
		particles.px += - kick_x * self.ds
		particles.py += - kick_y * self.ds

		## recompute relativistic momentum argument with new transverse momenta
		argumnent = 1. / self.compute_momentum_arg(particles)

		## second half drift
		particles.x += - particles.px * relativistic_factor * argumnent * self.ds / 2.
		particles.y += - particles.py * relativistic_factor * argumnent * self.ds / 2.


		return


	


