

import numpy as np
import matplotlib.pyplot as plt

import colormaps as cmaps

import solvers
import particles
import analytic_models

pi = np.pi
sigma_x = 0.005
x_max = 0.05
Q = 1.0e-9
x = np.linspace(-x_max,x_max,100)
y = 0 * x

gaussian = analytic_models.uniform_gaussian(sigma_r = sigma_x, Q_0 = Q)
r,phi = gaussian.compute_phi(x,y)

#kv = analytic_models.kv(r_0 = sigma_x, Q_0 = Q)
#r,phi = kv.compute_phi(x,y)



## Solver parameters
lambda_0 = 20.0 * x_max
n = 40

solver = solvers.field_solver_2D(lambda_x_0 = lambda_0, lambda_y_0 = lambda_0, 
	n_modes_x = n, n_modes_y = n )

## Particle distributions
particles = particles.particles_2D(Q0 = Q, N = 50000)
particles.construct_gaussian(sigma_x = sigma_x, sigma_y = sigma_x)
#particles.construct_kv(r_0 = sigma_x)


## Solver actions
solver.compute_phi(particles)
phi_grid, XX, YY = solver.compute_phi_mesh(xmax = x_max, ymax = x_max)


## post processing and  plotting 
phi_num = phi_grid[:,50] / np.max(phi_grid[:,50]) # - phi_grid[0,50])
phi_a = phi / np.max(phi) # - phi[0])

plt.figure()
plt.hexbin(particles.x,particles.y,cmap = cmaps.viridis) #plt.cm.RdBu)

plt.figure()
plt.plot(YY[:,50],phi_num, label = 'spectral solver')
plt.plot(x,phi_a, label = 'analytic model')
plt.legend(loc = 0)

plt.show()
