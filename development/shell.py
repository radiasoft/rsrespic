

import numpy as np
import matplotlib.pyplot as plt

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
n_particles = 100

## Field solver parameters 
L_0 = 15. * sigma_x ## Half the domain size
L_min = L_0 / 100 ## minimum wavelength to resolve

## Grid definitions for plotting 
x_max = 1.2 * sigma_x ## Half the plotting region
n_grid = 41
grid_index = int(n_grid / 2) 

## Here we construct the analytic solution for comparison purposes 
x = np.linspace(-x_max,x_max,50)
y = 0 * x
gaussian = analytic_models.uniform_gaussian(sigma_r = sigma_x, Q_0 = Q)
r_g,phi_g = gaussian.compute_phi(x,y)

## This is where we initialize a gaussian distribuiton
distribution = particles.distribution(N = n_particles)
distribution.construct_kv(r_0 = sigma_x)

## Particle distributions
## I am consturcting both tent and delta functions for comparison purposes 
particles_tent = particles.particles_2D_tent(dx_tent = 2. * L_min, dy_tent = 2. * L_min 
	, Q_0 = Q, N = n_particles)

particles_delta = particles.particles_2D_delta(Q_0 = Q, N = n_particles)

## Here the particle distributions are initialized 
particles_tent.initialize_particles(distribution)
particles_delta.initialize_particles(distribution)


## Define the fields 
fields_tent = fields.cartesian_2D(L_x = L_0, L_y = L_0,
	L_x_min = L_min, L_y_min = L_min)
fields_delta = fields.cartesian_2D(L_x = L_0, L_y = L_0,
	L_x_min = L_min, L_y_min = L_min)

## This is where we instantiate the solver
field_solver = solvers.field_solver_2D()

## Compute the potentials at the particle coordinates
field_solver.compute_phi(fields_tent, particles_tent)
field_solver.compute_phi(fields_delta, particles_delta)

## Compute the potentials on a grid
field_solver.compute_phi_mesh(fields_tent, xmax = x_max, ymax = x_max, n_grid = n_grid)
field_solver.compute_phi_mesh(fields_delta, xmax = x_max, ymax = x_max, n_grid = n_grid)


## post processing and  plotting 

plot = True

if plot:

	plt.figure()

	plt.subplot(1,2,1)
	plt.plot(fields_tent.y_grid[:, grid_index] / sigma_x, 
		fields_tent.phi_grid[:, grid_index], label = 'tent functions', linewidth = 2.)

	plt.plot(fields_delta.y_grid[:, grid_index] / sigma_x, 
		fields_delta.phi_grid[:, grid_index], label = 'delta functions', linewidth = 2.)

	#plt.plot(x / sigma_x, phi_g, label = 'analytic', linewidth = 2.)

	plt.legend(loc = 0)
	plt.xlabel(r'Position [$\sigma_x$]')
	plt.ylabel('Potential [V]')

	plt.subplot(1,2,2)
	plt.plot(fields_tent.y_grid[:, grid_index] / sigma_x,
		fields_tent.phi_grid[:, grid_index] - fields_delta.phi_grid[:, grid_index], linewidth = 2.)
	plt.xlabel(r'Position [$\sigma_x$]')
	plt.ylabel(r'$\Phi_\mathrm{tent} - \Phi_\mathrm{\detla}$ [V]')


	plt.show()
