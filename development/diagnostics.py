
import numpy as np
import matplotlib.pyplot as plt


class bunch_statistics:

    def __init__(self, divergence_coordinates = False):

        self.s = []

        self.step_data = []

        self.plot_keys = {'s': 0, 'rms_x': 1, 'rms_xp': 2, 
            'ex_rms': 3, 'rms_y': 4, 'rms_yp': 5, 'ey_rms': 6, 'r_beam': 7,
            'beta_x': 8, 'alpha_x': 9, 'beta_y': 10, 'alpha_y': 11 }


        self.divergence_coordinates = divergence_coordinates


    def update(self, s, particles):

        self.s.append(s)

        x = particles.x / 100.
        y = particles.y / 100.

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

        r_beam = np.max(np.sqrt(particles.x ** 2 + particles.y **2))

        ex = np.sqrt(np.dot(x, x) * np.dot(xp, xp) - np.dot(xp, x)**2) / len(xp)
        ey = np.sqrt(np.dot(y, y) * np.dot(yp, yp) - np.dot(yp, y)**2) / len(yp)

        beta_x = np.dot(x, x) / len(x) / ex
        beta_y = np.dot(y, y) / len(y) / ey
        alpha_x = -np.dot(xp, x) / len(x) / ex 
        alpha_y = -np.dot(yp, y) / len(y) / ey 


        step_data = [s, x_rms, xp_rms, ex, y_rms, yp_rms, ey, r_beam,
        beta_x, alpha_x, beta_y, alpha_y]


        self.step_data.append(step_data)


    def get_parameter(self, key):
        parameters = np.asarray(self.step_data)

        return parameters[:, self.plot_keys[key]]


    def save(self, file_name):
        all_data = np.asarray(self.step_data)
        np.save(file_name, all_data)

        return

    def load(self, file_name):
        all_data = np.load(file_name)
        self.step_data = all_data

        return

    def plot(self, x_key, y_key):
        all_data = np.asarray(self.step_data)

        x = all_data[:,self.plot_keys[x_key]]
        y = all_data[:,self.plot_keys[y_key]]

        plt.plot(x, y, linewidth = 2.)


