import numpy as np

from constants import cgs_constants
from numpy import exp, sin, einsum

pi = np.pi
q = cgs_constants['q']
c = cgs_constants['c'] 


## Convert units to cgs from mks

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

		trash, kx_mat = np.meshgrid(particles.x, fields.k_x_vector) ## 1/cm
		trash, ky_mat = np.meshgrid(particles.y, fields.k_y_vector) ## 1/cm 

		exp_x = np.exp(1j * kx1) * particles.lambda_twiddle(kx_mat, particles.x_extent) / np.sqrt(fields.lambda_x_0) ## no units
		exp_y = np.exp(1j * ky1) * particles.lambda_twiddle(ky_mat, particles.y_extent) / np.sqrt(fields.lambda_y_0) ## no units
	
		ptcl_exponential = np.einsum('mp, np -> mn', exp_x, exp_y) ## no units

		phi = einsum('xy, xy -> xy', ptcl_exponential, fields.k_sq_inv) ## no units

		fields.phi = - phi * particles.charge * particles.weight  * 4. * pi / ((particles.gamma ** 2 * particles.beta) ) ## statC s / cm
		
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

		phi = fields.phi 
		#statcolomb s / cm

		kx4 = np.einsum('m,i -> mi', fields.k_x_vector, xarray) ## no units
		ky4 = np.einsum('n,j -> nj', fields.k_y_vector, yarray) ## no units

		exp_x = np.exp(-1j * kx4) ## no units
		exp_y = np.exp(-1j * ky4) ## no units

		#Field components are exp(-i(kxx + kyy))
		phi_modes = np.einsum('mi, nj -> mnij', exp_x, exp_y) ## no units

		#now multiply by the sigma matrix, component-wise - using shared mn indices
		phi_vals = np.einsum('mn,mnij->ij',phi, phi_modes) ## statC s / cm

		# statcolomb s / cm
		fields.phi_grid = phi_vals - np.min(phi_vals)
		fields.x_grid = XX
		fields.y_grid = YY

		return

	def compute_psi_particles(self, fields, particles):

		phi = fields.phi
		## statC s / cm

		kx4 = np.einsum('m,i -> mi', fields.k_x_vector, particles.x) ## no units
		ky4 = np.einsum('n,i -> ni', fields.k_y_vector, particles.y) ## no units

		exp_x = np.exp(-1j * kx4) ## no units
		exp_y = np.exp(-1j * ky4) ## no units

		modes = np.einsum('mp, np -> mnp', exp_x, exp_y) ## no units

		psi_vals = np.einsum('mn, mnp -> p', phi, modes)

		fields.psi_vals = psi_vals

		return



	def compute_kick(self, fields, particles):

		phi = fields.phi
		## statC s / cm

		kx4 = np.einsum('m,i -> mi', fields.k_x_vector, particles.x) ## no units
		ky4 = np.einsum('n,i -> ni', fields.k_y_vector, particles.y) ## no units

		exp_x = np.exp(-1j * kx4) ## no units
		exp_y = np.exp(-1j * ky4) ## no units

		kick_modes = np.einsum('mp, np -> mnp', exp_x, exp_y) ## no units

		grad_psi_x = np.einsum('mn, m, mnp -> p', phi, -1j * fields.k_x_vector, kick_modes) ## statC s / cm^2
		grad_psi_y = np.einsum('mn, n, mnp -> p', phi, -1j * fields.k_y_vector, kick_modes) ## statC s / cm^2

		fields.kick_x = np.real(grad_psi_x) * particles.charge * particles.weight / (particles.beta * c) ## statC^2 s^2 / cm^3  
		fields.kick_y = np.real(grad_psi_y) * particles.charge * particles.weight / (particles.beta * c) ## statC^2 s^2 / cm^3  

		return


class kinetics_solver_SC2D:

	def __init__(self, ds = 1.0e-2):

		self.ds = ds



	def compute_momentum_arg(self,particles):

		## Calculate the square root term in the hamaltonian
		argument = np.sqrt(particles.beta**2 * particles.p_xi **2 - particles.px **2 - particles.py**2 - (particles.m_0 * particles.weight)**2 * c**2)
		
		return argument


	def step(self, particles, fields, field_solver):

		## The first thing is to trasform the longitudional momentum into canaonical coordiantes
		particles.compute_p_xi()

		## This computes the square root term of the hamaltonian that is used in all drifts
		argument = 1. / self.compute_momentum_arg(particles)
		relativistic_factor = 1. / particles.beta ** 2 - 1. 

		## First half drift
		particles.x += - particles.px * relativistic_factor * argument * self.ds / 2. 
		particles.y += - particles.py * relativistic_factor * argument * self.ds / 2. 

		## compute psi
		field_solver.compute_phi(fields, particles) 

		## compute the kick
		field_solver.compute_kick(fields, particles) ## statC s / cm^2

		## get the kick
		kick_x = fields.kick_x  ## statC^2 s^2 / cm^3 
		kick_y = fields.kick_y  ## statC^2 s^2 / cm^3 

		## apply the kick 
		particles.px += - kick_x * self.ds ##statC^2 s^2 / cm^2 --> This should be statC^2 s / cm^2
		particles.py += - kick_y * self.ds ##statC^2 s^2 / cm^2

		## recompute relativistic momentum argument with new transverse momenta
		argument = 1. / self.compute_momentum_arg(particles)

		## second half drift
		particles.x += - particles.px * relativistic_factor * argument * self.ds / 2.
		particles.y += - particles.py * relativistic_factor * argument * self.ds / 2.


		return


	


class kinetics_solver_FODO_Cell:

	def __init__(self, ds = 1.0e-2):

		self.ds = ds



	def compute_momentum_arg(self,particles):

		## Calculate the square root term in the hamaltonian
		argumnent = np.sqrt(particles.beta**2 * particles.p_xi **2 - particles.px **2 - particles.py**2 - particles.m_0**2 * c**2)
		
		return argumnent


	def step(self, particles, fields, field_solver):



		return

