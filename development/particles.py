import numpy as np
import constants
from constants import cgs_constants
from scipy.special import erf

pi = np.pi

## Convert units to cgs from mks
elementary_charge = cgs_constants['q']
c = cgs_constants['c']
m_e = cgs_constants['m_e']

class distribution:

	def __init__(self, N = 0):

		self.N = N

	
	def construct_uniform_guassian_2D(self, x0 = 0, y0 =0, xp0 = 0, yp0 = 0,
			sigma_x = 0, sigma_y = 0, sigma_xp = 0, sigma_yp = 0):
		

		## note xp amd yp are in radians. 

		sigma_x = sigma_x
		sigma_y = sigma_y
		sigma_xp = sigma_xp
		sigma_yp = sigma_yp

		sigma_xxp = 0.0
		sigma_xy = 0.0
		sigma_xyp = 0.0
		sigma_xpy = 0.0
		sigma_xpyp = 0.0
		sigma_yyp = 0.0
			
		mean = [x0, xp0, y0, yp0]

		cov = [[sigma_x**2, sigma_xxp**2 , sigma_xy**2 ,sigma_xyp**2],
			[sigma_xxp**2,sigma_xp**2, sigma_xpy**2 ,sigma_xpyp**2],
			[sigma_xy**2,sigma_xpy**2,sigma_y**2,sigma_yyp**2],
			[sigma_xyp**2,sigma_xpyp**2,sigma_yyp**2,sigma_yp**2]]

		x, xp, y, yp = np.random.multivariate_normal(mean, cov, self.N).T
		
		z = np.zeros(len(x))
		pz = np.zeros(len(x))
		self.x = x
		self.y = y
		self.xp = xp
		self.yp = yp
		self.z = z
		self.pz = pz

		self.type = 'gaussian'
		return


	def construct_kv(self, x_a = 0.1 , x_b = 0.1, y_a = 0.1, y_b = 0.1):
		
		self.e_x = np.pi * x_a * x_b / 100.
		self.e_y = np.pi * y_a * y_b / 100.

		r_x = np.random.uniform(0, x_a**2, self.N)
		theta_x = np.random.uniform(0, 2. * np.pi, self.N)

		r_y = np.random.uniform(0, y_a**2, self.N)
		theta_y = np.random.uniform(0, 2. * np.pi, self.N) 

		x = np.sqrt(r_x) * np.cos(theta_x)
		xp = np.sqrt(r_x) * np.sin(theta_x) * x_b / x_a

		y = np.sqrt(r_y) * np.sin(theta_y)
		yp = np.sqrt(r_y) * np.cos(theta_y) * y_b / y_a

		#px = np.zeros(len(x))
		#py = np.zeros(len(x))
		self.type = 'KV'

		z = np.zeros(len(x))
		pz = np.zeros(len(x))

		self.x = x
		self.y = y
		self.xp = xp
		self.yp = yp
		self.z = z
		self.pz = pz

		return



class particles_2D_delta:

	def __init__(self, distribution, bunch_charge = 1.0, 
		species_charge = elementary_charge, species_mass = m_e, K_e = 1.0e6):	
		
		self.x_extent = 1.
		self.y_extent = 1.

		self.bunch_charge = bunch_charge
		self.charge = species_charge
		self.m_0 = species_mass

		self.mc2 = (self.m_0 / 1000.) * (c / 100.) ** 2 / (constants.charge_cgs_to_mks(self.charge))

		self.gamma = K_e / self.mc2 + 1.

		self.weight = bunch_charge / species_charge / distribution.N

		self.N = distribution.N

		self.set_gamma(self.gamma)
		self.initialize_particles(distribution)


	def set_gamma(self,gamma):

		self.gamma = gamma
		self.beta = np.sqrt( 1. - 1. / gamma **2 )

	def initialize_particles(self,distribution):

		self.pz = self.gamma * self.m_0 * c * self.weight * self.beta
		self.pt = distribution.pz + self.gamma * self.m_0 * c * self.weight
		self.x = distribution.x
		self.px = distribution.xp * self.pz
		self.y = distribution.y
		self.py = distribution.yp * self.pz
		self.z = distribution.z

		if distribution.type == 'KV':
			self.e_x = distribution.e_x
			self.e_y = distribution.e_y
		else:
			self.e_x = np.sqrt(np.dot(distribution.x, distribution.x) * np.dot(distribution.xp, distribution.xp) - np.dot(distribution.x, distribution.xp)) / distribution.N
			self.e_y = np.sqrt(np.dot(distribution.y, distribution.y) * np.dot(distribution.yp, distribution.yp) - np.dot(distribution.y, distribution.yp)) / distribution.N

		self.compute_p_xi()


	def compute_p_xi(self):

		self.p_xi = self.pt / self.beta

		return


	def lambda_twiddle(self, k_matrix, extent):
		
		row,col = k_matrix.shape

		return np.ones(k_matrix.shape)


	def write_opal_distribution(self, z_extent = 1.0):
		## This will write a file that can be loaded into opal with an extent of 1m 

		## Convert units on particles
		x_out = self.x / 100.
		y_out = self.y / 100.
		z_out = np.random.uniform(0, z_extent, self.N)
		px_out = self.px / (self.m_0 * c * self.weight)
		py_out = self.py / (self.m_0 * c * self.weight)
		pz_out = self.pz / (self.m_0 * c * self.weight)

		particle_array = np.asarray([x_out, px_out, y_out, py_out, z_out, pz_out])

		np.savetxt('opal_input_distribution.txt', particle_array, delimiter = '\t')

		return

class particles_2D_tent:

	def __init__(self, distribution, dx_tent = 0.1, dy_tent = 0.1,  bunch_charge = 1.0, 
		species_charge = elementary_charge, species_mass = m_e, K_e = 1.0e6):

		self.x_extent = dx_tent
		self.y_extent = dy_tent
		
		self.bunch_chage = bunch_charge
		self.charge = species_charge
		self.m_0 = species_mass

		self.mc2 = (self.m_0 / 1000.) * (c / 100.) ** 2 / (constants.charge_cgs_to_mks(self.charge))

		self.gamma = K_e / self.mc2 + 1

		self.weight = bunch_charge / species_charge / distribution.N

		self.N = distribution.N

		self.set_gamma(self.gamma)
		self.initialize_particles(distribution)


	def set_gamma(self,gamma):

		self.gamma = gamma
		self.beta = np.sqrt( 1. - 1. / gamma **2 )


	def initialize_particles(self,distribution):

		self.pz = self.gamma * self.m_0 * c * self.weight * self.beta
		self.pt = distribution.pz + self.gamma * self.m_0 * c * self.weight
		self.x = distribution.x
		self.px = distribution.xp * self.pz
		self.y = distribution.y
		self.py = distribution.yp * self.pz
		self.z = distribution.z

		if distribution.type == 'KV':
			self.e_x = distribution.e_x
			self.e_y = distribution.e_y
		else:
			self.e_x = np.sqrt(np.dot(distribution.x, distribution.x) * np.dot(distribution.xp, distribution.xp) - np.dot(distribution.x, distribution.xp)) / distribution.N
			self.e_y = np.sqrt(np.dot(distribution.y, distribution.y) * np.dot(distribution.yp, distribution.yp) - np.dot(distribution.y, distribution.yp)) / distribution.N

		self.compute_p_xi()


	def compute_p_xi(self):

		self.p_xi = self.pt / self.beta

		return


	def lambda_twiddle(self, k, tent_width):

		arg = k * tent_width / (2. * pi) 
 
		lambda_twiddle = np.sinc(arg) ** 2

		return lambda_twiddle



