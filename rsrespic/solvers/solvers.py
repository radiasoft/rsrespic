import numpy as np

from rsrespic.utilities import constants 
from numpy import exp, sin, einsum

pi = np.pi
q = constants.cgs_constants['q']
c = constants.cgs_constants['c'] 


## Convert units to cgs from mks


class sine_transform_2D(object):

	def __init__(self):

		self.name = '2-d electrostatic solver using sin transform'


	def compute_grad_psi(self, fields, particles):


		kappa_1 = particles.bunch_charge / particles.N

		M = fields.n_modes_x
		N = fields.n_modes_y
		Lx = fields.L_x
		Ly = fields.L_y
		
		## shift the particles into the domain for computing the fiels
		x = particles.x + Lx / 2.
		y = particles.y + Ly / 2.

		arr_x = np.exp(1.j *np.pi *x /Lx) 
		arr_y = np.exp(1.j *np.pi *y /Ly) 

		E_x = 0.0 *x  # field components at the locations of the particles 
		E_y = 0.0 *y 

		for i_m in np.arange(M): 
			m = i_m +1 
			tmp_x = arr_x**m 
			tx_r = np.real(tmp_x) 
			tx_i = np.imag(tmp_x) 
			for i_n in np.arange(N): 
				n = i_n +1 
				tmp_y = arr_y**n 
				ty_r = np.real(tmp_y) 
				ty_i = np.imag(tmp_y) 
				c_mn_rho = np.sum(tx_i  *ty_i) 
				c_mn_rho *= 4./ (Lx *Ly) 
				c_mn_phi = 4. *np.pi * kappa_1 *c_mn_rho / ( (m*np.pi/Lx)*(m*np.pi/Lx) +(n*np.pi/Ly)*(n*np.pi/Ly) )  # "-" from the Laplacian, "-" from the -4pi in the rhs 
				E_x -= (m*np.pi/Lx) *c_mn_phi *tx_r *ty_i  #  *np.cos(m *np.pi *x[ix] /Lx) *np.sin(n *np.pi *y[iy] /Ly)  # -d/dx 
				E_y -= (n*np.pi/Ly) *c_mn_phi *tx_i *ty_r  #  *np.sin(m *np.pi *x[ix] /Lx) *np.cos(n *np.pi *y[iy] /Ly)  # -d/dy
		
		fields.psi_x = E_x / (particles.gamma ** 2.)  ## statC^2 s^2 / cm^3  
		fields.psi_y = E_y / (particles.gamma ** 2.) ## statC^2 s^2 / cm^3  


class field_solver_2D(object):

	def __init__(self):

		self.name = '2-d electrostatic field solver'

	def compute_mode_coefficients(self, fields, particles):

		## Setup the coefficients for the 
		kx1 = np.einsum('m, p -> mp', fields.k_x_vector, particles.x)
		ky1 = np.einsum('n, p -> np', fields.k_y_vector, particles.y)

		trash, kx_mat = np.meshgrid(particles.x, fields.k_x_vector) ## 1/cm
		trash, ky_mat = np.meshgrid(particles.y, fields.k_y_vector) ## 1/cm 

		exp_x = np.exp(1j * kx1) * particles.lambda_twiddle(kx_mat, particles.x_extent) / fields.lambda_x_0 ## no units
		exp_y = np.exp(1j * ky1) * particles.lambda_twiddle(ky_mat, particles.y_extent) / fields.lambda_y_0 ## no units

		ptcl_exponential = np.einsum('mp, np -> mn', exp_x, exp_y) ## no units

		unscalled_coefficients = einsum('xy, xy -> xy', ptcl_exponential, fields.k_sq_inv) ## no units

		fields.mode_coefficients = - unscalled_coefficients * particles.weight  * 4. * pi  * q * np.sqrt(2)/ ((particles.gamma ** 2) ) ## statC s / cm
		
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

		phi = fields.mode_coefficients
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

		phi = fields.mode_coefficients 
		## statC s / cm

		kx4 = np.einsum('m,i -> mi', fields.k_x_vector, particles.x) ## no units
		ky4 = np.einsum('n,i -> ni', fields.k_y_vector, particles.y) ## no units

		exp_x = np.exp(-1j * kx4) ## no units
		exp_y = np.exp(-1j * ky4) ## no units

		modes = np.einsum('mp, np -> mnp', exp_x, exp_y) ## no units

		psi_vals = np.einsum('mn, mnp -> p', phi, modes)

		fields.psi_vals = psi_vals

		return



	def compute_grad_psi(self, fields, particles):

		phi = fields.mode_coefficients 
		## statC s / cm

		kx4 = np.einsum('m,i -> mi', fields.k_x_vector, particles.x) ## no units
		ky4 = np.einsum('n,i -> ni', fields.k_y_vector, particles.y) ## no units

		exp_x = np.exp(-1j * kx4) ## no units
		exp_y = np.exp(-1j * ky4) ## no units

		kick_modes = np.einsum('mp, np -> mnp', exp_x, exp_y) ## no units

		grad_psi_x = np.einsum('mn, m, mnp -> p', phi, -1j * fields.k_x_vector, kick_modes) ## statC s / cm^2
		grad_psi_y = np.einsum('mn, n, mnp -> p', phi, -1j * fields.k_y_vector, kick_modes) ## statC s / cm^2

		fields.psi_x = np.real(grad_psi_x) ## statC^2 s^2 / cm^3  
		fields.psi_y = np.real(grad_psi_y)  ## statC^2 s^2 / cm^3  

		return




class symplectic_maps:

	def __init__(self, name = 'maps'):

		self.name = name 

	def space_charge_kick_2D_sine(self, fields, particles, ds = 0.0):

		## compute the fields from the particles using the sin transform
		fields.solver.compute_grad_psi(fields, particles)

		## compute the kick 
		kick_x = fields.psi_x * particles.charge  * particles.weight / (particles.beta * c)  ## statC^2 s^2 / cm^3 
		kick_y = fields.psi_y * particles.charge  * particles.weight / (particles.beta * c)  ## statC^2 s^2 / cm^3 

		## apply the kick 
		particles.px = particles.px + kick_x * ds ##statC^2 s^2 / cm^2 --> This should be statC^2 s / cm^2
		particles.py = particles.py + kick_y * ds ##statC^2 s^2 / cm^2
		particles.pt = particles.pt +  0. 



	def space_charge_kick_2D(self, fields, particles, ds = 0.0):
		
		## compute mode coefficients
		fields.solver.compute_mode_coefficients(fields, particles) 

		## compute the field gradients
		fields.solver.compute_grad_psi(fields, particles) ## statC s / cm^2

		## compute the kick 
		kick_x = fields.psi_x * particles.charge  * particles.weight / (particles.beta * c)  ## statC^2 s^2 / cm^3 
		kick_y = fields.psi_y * particles.charge  * particles.weight / (particles.beta * c)  ## statC^2 s^2 / cm^3 

		## apply the kick 
		particles.px = particles.px + kick_x * ds ##statC^2 s^2 / cm^2 --> This should be statC^2 s / cm^2
		particles.py = particles.py + kick_y * ds ##statC^2 s^2 / cm^2
		particles.pt = particles.pt +  0. 

		#return fields, particles


	def drift(self, particles, ds = 0.0):
		## Compute the normalizing factor for the momentum
		argument = np.sqrt( (particles.beta * particles.p_xi) **2 - particles.px **2 
			- particles.py**2 - (particles.m_0 * particles.weight * c)**2)

		## compute the rate of drift per s 
		dx_ds = particles.px / argument 
		dy_ds = particles.py / argument
		dz_ds = 0.

		## update particle positoins 
		particles.x = particles.x + particles.px / argument * ds
		particles.y = particles.y + particles.py / argument * ds
		particles.z = particles.z + 0



	def thin_quad(self, fields, particles, ds = 0.0, kappa = 0.0):

		## Thin quad has no drift 
		dx_ds = 0.
		dy_ds = 0.
		dz_ds = 0.

		## assume a positron convention here kappa is the same as elegant K1
		dpx = kappa * particles.x * ds * (particles.pz / 10000.)
		dpy = - kappa * particles.y * ds * (particles.pz / 10000.)
		dpz = 0.

		particles.px += dpx
		particles.py += dpy
		particles.pt += 0. 


	def reset_modes(self, fields, particles):

		sigma_x = np.std(particles.x)
		sigma_y = np.std(particles.y)
		#fields.L_y = 

		return


