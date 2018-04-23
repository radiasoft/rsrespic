import numpy as np
from scipy.constants import c, mu_0, m_e
from scipy.constants import epsilon_0 as e_0
from scipy.constants import elementary_charge as q
from scipy.special import erf


pi = np.pi



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
		
		self.x = x
		self.y = y
		self.px = px
		self.py = py

		return x,px,y,py


	def construct_kv(self, r_0 = 1.0):
		
		r = np.random.uniform(0, r_0**2, self.N)
		theta = np.random.uniform(0, 2*pi, self.N)

		x = np.sqrt(r) * np.cos(theta)
		y = np.sqrt(r) * np.sin(theta)

		px = np.zeros(len(x))
		py = np.zeros(len(x))
		
		self.x = x
		self.y = y
		self.px = px
		self.py = py
		
		return x,px,y,py



class particles_2D_delta:

	def __init__(self, Q_0 = 1.0, N = 1000, gamma = 1.0, m = 1.):	

		self.Q_0 = Q_0
		self.N = N
		self.w = self.Q_0 / N 
		self.x = np.zeros(N)
		self.px = np.zeros(N)
		self.y = np.zeros(N)
		self.py = np.zeros(N)
		self.x_extent = 1.
		self.y_extent = 1.
		self.gamma = gamma
		self.beta = np.sqrt( 1. - 1. / gamma **2 )
		self.m = m
		self.pz = 10.

	def initialize_particles(self,distribution):

		self.x = distribution.x
		self.px = distribution.px
		self.y = distribution.y
		self.py = distribution.py


	def lambda_twiddle(self, k_matrix, extent):
		
		row,col = k_matrix.shape

		return np.ones(k_matrix.shape)




class particles_2D_tent:

	def __init__(self, dx_tent = 0.1, dy_tent = 0.1, Q_0 = 1.0, N = 1000, gamma = 1.0, m = 1.):
		#self.coordinates = coordinates
		#self.n_particles = len(coordinates)	
		self.x_extent = dx_tent
		self.y_extent = dy_tent
		self.Q_0 = Q_0
		self.N = N
		self.w = self.Q_0 / N
		self.x = np.zeros(N)
		self.px = np.zeros(N)
		self.y = np.zeros(N)
		self.py = np.zeros(N)
		self.gamma = gamma
		self.beta = np.sqrt( 1. - 1. / gamma **2 )
		self.m = m
		self.pz = 10.

	def initialize_particles(self,distribution):

		self.x = distribution.x
		self.px = distribution.px
		self.y = distribution.y
		self.py = distribution.py


	def lambda_twiddle(self, k, tent_width):

		arg = k * tent_width / (2. * pi) 
 
		lambda_twiddle = np.sinc(arg) ** 2

		return lambda_twiddle



