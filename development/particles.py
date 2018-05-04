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

		x, px, y, py = np.random.multivariate_normal(mean, cov, self.N).T
		
		z = np.zeros(len(x))
		pz = np.zeros(len(x))

		self.x = x
		self.y = y
		self.px = px
		self.py = py
		self.z = z
		self.pz = pz

		return


	def construct_kv(self, r_0 = 1.0, x0 = 0.0, y0 = 0.0, pr_0 = 0, pth_0 = 0):
		
		r = np.random.uniform(0, r_0**2, self.N)
		theta = np.random.uniform(0, 2*pi, self.N)

		x = np.sqrt(r) * np.cos(theta) + x0
		y = np.sqrt(r) * np.sin(theta) + y0

		px = np.zeros(len(x))
		py = np.zeros(len(x))
		
		z = np.zeros(len(x))
		pz = np.zeros(len(x))

		self.x = x
		self.y = y
		self.px = px
		self.py = py
		self.z = z
		self.pz = pz

		return



class particles_2D_delta:

	def __init__(self, distribution, bunch_charge = 1.0, 
		species_charge = elementary_charge, species_mass = m_e, K_e = 1.0e6):	
		
		self.x_extent = 1.
		self.y_extent = 1.

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

		self.x = distribution.x
		self.px = distribution.px
		self.y = distribution.y
		self.py = distribution.py
		self.z = distribution.z
		self.pt = distribution.pz + self.gamma * self.m_0 * c * self.weight

		self.compute_p_xi()


	def compute_p_xi(self):

		self.p_xi = self.pt / self.beta

		return


	def lambda_twiddle(self, k_matrix, extent):
		
		row,col = k_matrix.shape

		return np.ones(k_matrix.shape)




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

		self.x = distribution.x
		self.px = distribution.px
		self.y = distribution.y
		self.py = distribution.py
		self.z = distribution.z
		self.pt = distribution.pz + self.gamma * self.m_0 * c * self.weight

		self.compute_p_xi()


	def compute_p_xi(self):

		self.p_xi = self.pt / self.beta

		return


	def lambda_twiddle(self, k, tent_width):

		arg = k * tent_width / (2. * pi) 
 
		lambda_twiddle = np.sinc(arg) ** 2

		return lambda_twiddle



