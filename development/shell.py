

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams.update({'font.size': 16})

import colormaps as cmaps

import solvers
import particles
import analytic_models

pi = np.pi
sigma_x = 0.05

L_0 = 20. * sigma_x ## Half the domain size
L_min = L_0 / 160. ## minimum wavelength to resolve

x_max = 1 * sigma_x ## Half the plotting region

Q = 50.0e-12
x = np.linspace(-x_max,x_max,100)
y = 0 * x

gaussian = analytic_models.uniform_gaussian(sigma_r = sigma_x, Q_0 = Q)
r_g,phi_g = gaussian.compute_phi(x,y)

kv = analytic_models.kv(r_0 = sigma_x, Q_0 = Q)
r_k,phi_k = kv.compute_phi(x,y)


## Solver parameters
solver = solvers.field_solver_2D(L_x = L_0, L_y = L_0,
	L_x_min = L_min, L_y_min = L_min)

## Particle distributions
particles = particles.particles_2D(Q0 = Q, N = 50000)
#particles.construct_gaussian(sigma_x = sigma_x, sigma_y = sigma_x)
particles.construct_kv(r_0 = sigma_x)

n_grid = 21
grid_index = int(n_grid / 2) 
## Solver actions
solver.compute_phi(particles)
phi_grid, XX, YY = solver.compute_phi_mesh(xmax = x_max, ymax = x_max, n_grid = n_grid)


## post processing and  plotting 
phi_num = phi_grid[:,grid_index] #/ np.max(phi_grid[:,50]) # - phi_grid[0,50])

plot = True


if plot:
	plt.figure()
	plt.plot(YY[:, grid_index]/ sigma_x, phi_num, label = 'spectral solver', linewidth = 2.)
	plt.plot(x / sigma_x, phi_k, label = 'analytic model kv', linewidth = 2.)
	plt.plot(x / sigma_x, phi_g, label = 'analytic model gaussian', linewidth = 2.)

	plt.legend(loc = 0)
	plt.xlabel(r'Position [$\sigma_x$]')
	plt.ylabel('Potential [V]')

	plt.figure()
	plt.contourf(XX / sigma_x, YY / sigma_x, phi_grid, 20, cmap = cmaps.viridis)
	plt.xlabel(r'Position [$\sigma_x$]')
	plt.ylabel(r'Position [$\sigma_y$]')
	plt.colorbar()

	plt.show()
