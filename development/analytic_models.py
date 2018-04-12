import numpy as np
from scipy.constants import c, mu_0, m_e
from scipy.constants import epsilon_0 as e_0
from scipy.constants import elementary_charge as q
from scipy.special import erf


pi = np.pi




class uniform_gaussian:
	""" Class for computing analytic solutions to fields and potentials 
	for a cylidrical gaussian besm""" 

	def __init__(self, sigma_r = 1, Q_0 = 1):

		""" Q_0 is total charge in C and sigma_r is the beam sigma in m"""
		self.sigma_r = sigma_r
		self.Q_0 = Q_0

	def compute_phi(self, x, y):	

		""" Convert from caetesian to cylindrical and compute potentials"""
		r = np.sqrt(x**2 + y**2)

		phi_r = self.Q_0 / (2 * pi * e_0 * r) * erf( r / (np.sqrt(2) * self.sigma_r) )

		return r, phi_r


	def compute_E(self, x, y):

		""" Convert from cartesian to cylindrtical and compute fields"""
		r = np.sqrt(x**2 + y**2)

		E_r = self.Q_0 / (e_0) * ( 1. / ( np.sqrt(2.) * pi**1.5 * self.sigma_r * r) * np.exp( - r**2 / ( 2. * self.sigma_r **2 )) -
			erf(r / (np.sqrt(2) * self.sigma_r)) / (2. * pi * r**2)) 

		return r,E_r



class waterbag:

	def __init__(self, r_0 = 0, Q_0 = 1):

		self.r_0 = r_0
		self.Q_0 = Q_0


	def compute_phi(self, x, y):

		
		
		return x



	def compute_E(self, x, y):

		return y
