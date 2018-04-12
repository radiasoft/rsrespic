

import numpy as np
import matplotlib.pyplot as plt

import colormaps as cmaps

from solvers import * 
from particles import *
from analytic_models import *

sigma_x = 0.1

gaussian = uniform_gaussian(sigma_r = sigma_x, Q_0 = 1.0)
#
x = np.linspace(-5.,5.,10000)
y = 0 * x

r,phi = gaussian.compute_phi(x,y)
#r,E_r = gaussian.compute_E(x,y)

#E = np.diff(phi) / np.mean(np.diff(r))



lambda_0 = 40.0
n = 40.0
solver = field_solver_2D(lambda_x_0 = lambda_0, lambda_y_0 = lambda_0, n_modes_x = n, n_modes_y = n )

particles = particles_2D()
particles.construct_gaussian_r(50000,sigma_x)
#particles.construct_kv(50000)
plt.figure()
plt.hexbin(particles.x,particles.y,cmap = cmaps.viridis) #plt.cm.RdBu)


solver.compute_phi(particles)

#plt.figure()
#plt.plot(solver.phi)

phi_grid, XX, YY = solver.compute_phi_mesh(xmax = 5.0, ymax = 5.0)

#print phi_grid

plt.figure()
plt.contourf(XX, YY, -phi_grid, 60 , cmap = cmaps.viridis) #plt.cm.RdBu)
plt.colorbar()

plt.figure()
plt.plot(YY[:,50],np.abs(phi_grid[:,50])/np.max(np.abs(phi_grid[:,50])), label = 'spectral solver')
plt.plot(x, np.exp(-x**2 / 2.))
plt.plot(x,phi/np.max(phi), label = 'analyticmodel')
plt.legend(loc = 0)

plt.show()
