

cgs_constants = {
'u': 1.660538782e10-24,
'mu_B':9.27400915e-21,
'a0':5.2917720859e-9,
'm_e':9.10938215e-28,
'q':4.80320427e-10,
'h_bar':1.0545716e10-27,
'c':2.99792458e10,
'k':1.3806504e-16,
'm_p':1.6726219e-24,
}



def charge_mks_to_cgs(q):

	return 2997919999.934 * q


def charge_cgs_to_mks(q):

	return q / 2997919999.934


def SV_to_V(SV):

	return SV*299.792458

def V_to_SV(V):

	return V/299.792458

def SVpercm_to_Vperm(E):

	return 29979.19999934*E

def Vperm_to_SVpercm(E):

	return E/29979.19999934