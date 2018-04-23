

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, mu_0, m_e

import time

import matplotlib
matplotlib.rcParams.update({'font.size': 16})

import colormaps as cmaps

import solvers
import particles
import fields
import analytic_models

pi = np.pi

## Particle definitions 
sigma_x = 0.05
Q = 50.0e-12
n_particles = 1000

## Field solver parameters 
L_0 = 15. * sigma_x ## Half the domain size
L_min = L_0 / 100 ## minimum wavelength to resolve
 


## This is where we initialize a gaussian distribuiton
distribution = particles.distribution(N = n_particles)
distribution.construct_kv(r_0 = sigma_x)

## Particle distributions
## I am consturcting both tent and delta functions for comparison purposes 
particles = particles.particles_2D_tent(dx_tent = 2. * L_min, dy_tent = 2. * L_min 
	, Q_0 = Q, N = n_particles, gamma = 10.0, m = m_e)

particles.initialize_particles(distribution)


## Define the fields 
fields = fields.cartesian_2D(L_x = L_0, L_y = L_0,
	L_x_min = L_min, L_y_min = L_min)

## This is where we instantiate the solver
field_solver = solvers.field_solver_2D()

# instantiate the kinetics solver
kinetics_solver = solvers.kinetics_solver_SC2D(ds = 1.0e-6)

r_max = []
s_beam = []

s = 0
k = 0
while k < 100:
	t0 = time.time()

	field_solver.compute_phi(fields, particles)
	field_solver.compute_kick(fields, particles)

	particles = kinetics_solver.step(particles, fields)

	r = np.sqrt(particles.x**2 + particles.y**2)
	r_max.append(np.max(r))

	print np.max(r)

	s = s + kinetics_solver.ds

	s_beam.append(s)
	t1 = time.time()
	k = k + 1
	print k, t1 - t0



plt.figure()
plt.plot(r_max)
plt.show()
