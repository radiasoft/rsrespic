import numpy as np
from scipy.constants import c, mu_0, m_e
from scipy.constants import epsilon_0 as e_0
from scipy.constants import elementary_charge as q
from scipy.special import erf
from scipy.special import gammainc
from scipy.special import gamma

pi = np.pi
gamma_em = 0.57721566490

def gamma_inc(x):

	return gamma(0.001)*(1 - gammainc(0.001,x))


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

		#phi_r = self.Q_0 / (2 * pi * e_0 * r) * erf( r / (np.sqrt(2) * self.sigma_r) )

		arg_r = r**2 / (2. * self.sigma_r ** 2)

		phi_r = (self.Q_0 * pi / (self.sigma_r ** 2 * e_0) * 1. / 2. * self.sigma_r**2 * (np.log( arg_r ) + gamma_inc(arg_r) + gamma_em))
		#phi_r = ( pi / (self.sigma_r ** 2 ) * 1. / 2. * self.sigma_r**2 * (np.log( arg_r ) + gamma_inc(arg_r) + gamma_em))

		return r, phi_r


	def compute_E(self, x, y):

		""" Convert from cartesian to cylindrtical and compute fields"""
		r = np.sqrt(x**2 + y**2)

		E_r = self.Q_0 / (e_0) * ( 1. / ( np.sqrt(2.) * pi**1.5 * self.sigma_r * r) * np.exp( - r**2 / ( 2. * self.sigma_r **2 )) -
			erf(r / (np.sqrt(2) * self.sigma_r)) / (2. * pi * r**2)) 

		return r,E_r

class kv:

	def __init__(self, r_0 = 1.0, Q_0 = 1.0):

		self.r_0 = r_0
		self.Q_0 = Q_0

	def compute_phi(self,x,y):

		r = np.sqrt(x**2 + y**2)

		phi_r = self.Q_0 * pi / (2. * e_0 * self.r_0**2) * ( r**2 * (r < self.r_0) + 
			(np.log(r) - np.log(self.r_0) + self.r_0 ** 2) * (r >= self.r_0) )

		return r,phi_r

	def compute_E(self,x,y):

		r = np.sqrt(x**2 + y**2)

		E_r = self.Q_0 / (e_0 * self.r_0**2) * ( r**2 * (r < self.r_0) + 
			(np.log(r) - np.log(self.r_0) + 1) * self.r_0**2 * (r >= self.r_0))



class waterbag:

	def __init__(self, r_0 = 0, Q_0 = 1):

		self.r_0 = r_0
		self.Q_0 = Q_0


	def compute_phi(self, x, y):

		
		
		return x



	def compute_E(self, x, y):

		return y
