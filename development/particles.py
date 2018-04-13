import numpy as np
from scipy.constants import c, mu_0, m_e
from scipy.constants import epsilon_0 as e_0
from scipy.constants import elementary_charge as q
from scipy.special import erf


pi = np.pi


	

class particles_2D:

	def __init__(self, dx_tent = 0.1, dy_tent = 0.1, Q0 = 1.0):
		#self.coordinates = coordinates
		#self.n_particles = len(coordinates)	
		self.dx_tent = dx_tent
		self.dy_tent = dy_tent
		self.Q0 = Q0

	def construct_gaussian(self, N, sigma_x = 1.0):
	
		w = self.Q0 / N

		sigma_x = sigma_x
		sigma_y = sigma_x
		sigma_xp = 0.01
		sigma_yp = 0.01


		sigma_xxp = 0.0005
		sigma_xy = 0.0000
		sigma_xyp = 0.0000
		sigma_xpy = 0.0000
		sigma_xpyp = 0.0000
		sigma_yyp = 0.0000
			
		mean = [0,0,0,0]


		cov = [[sigma_x**2, sigma_xxp**2 , sigma_xy**2 ,sigma_xyp**2],
			[sigma_xxp**2,sigma_xp**2, sigma_xpy**2 ,sigma_xpyp**2],
			[sigma_xy**2,sigma_xpy**2,sigma_y**2,sigma_yyp**2],
			[sigma_xyp**2,sigma_xpyp**2,sigma_yyp**2,sigma_yp**2]]

		x, px, y, py = np.random.multivariate_normal(mean, cov, N).T

		self.x = x
		self.y = y
		self.charge = w * np.ones(len(x))

	def construct_gaussian_r(self, N, sigma_r):

		r = np.abs(np.random.normal(0, sigma_r**2, N))
		theta = np.random.uniform(0,2*pi, N)

		x = np.sqrt(r) * np.cos(theta)
		y = np.sqrt(r) * np.sin(theta)

		self.x = x
		self.y = y
		self.charge = q * np.ones(len(x))

	def lambda_twiddle(self, k_x_vector, k_y_vector):

		return


	def construct_kv(self,N, r_0 = 0.5):

		r = np.random.uniform(0,r_0**2,N)
		theta = np.random.uniform(0,2*pi, N)

		x = np.sqrt(r) * np.cos(theta)
		y = np.sqrt(r) * np.sin(theta)

		self.x = x
		self.y = y
		self.charge = np.ones(len(x))


