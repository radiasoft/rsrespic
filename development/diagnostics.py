
import numpy as np
import matplotlib.pyplot as plt


class bunch_statistics:

    def __init__(self, divergence_coordinates = False):

        self.s = []

        self.step_data = []

        self.plot_keys = {'s': 0, 'rms_x': 1, 'rms_xp': 2, 
            'ex_rms': 3, 'rms_y': 4, 'rms_yp': 5, 'ey_rms': 6}

        self.divergence_coordinates = divergence_coordinates


    def update(self, s, particles):

        self.s.append(s)

        x_rms = np.std(particles.x)
        y_rms = np.std(particles.y)

        if self.divergence_coordinates:
            xp = particles.xp
            yp = particles.yp
        
        else:
            xp = particles.px / particles.pz
            yp = particles.py / particles.pz
        

        xp_rms = np.std(xp)
        yp_rms = np.std(yp)

        ex = np.sqrt(np.dot(particles.x,particles.x) * np.dot(xp,xp) - np.dot(xp, particles.x)**2) / len(xp)
        ey = np.sqrt(np.dot(particles.y,particles.y) * np.dot(yp,yp) - np.dot(yp, particles.y)**2) / len(yp)

        step_data = [s, x_rms, xp_rms, ex, y_rms, yp_rms, ey]

        self.step_data.append(step_data)


    def get_parameters(self, s):

        return s



    def plot(self, x_key, y_key):
        all_data = np.asarray(self.step_data)

        x = all_data[:,self.plot_keys[x_key]]
        y = all_data[:,self.plot_keys[y_key]]

        plt.plot(x, y, linewidth = 2.)


