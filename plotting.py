#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import math


def load_text_file(name):
    #This Function loads and plots specified result data
    print("Loading text file", name)
    #Read in data and sort it to get the Converged data file
    times = os.listdir("postProcessing/sets")
    times.sort()
    latest = times[-1]
    #Select the path to the specific file
    filepath = os.path.join("postProcessing", "sets", latest, 
                            name)
    filepath2 = os.path.join("postProcessing", "sets", latest, 
                           "hub_epsilon_k_p.xy")
    filepath3 = os.path.join("postProcessing", "sets", latest, 
                            "data_1D_U.xy")
    filepath4 = os.path.join("postProcessing", "sets", latest, 
                            "data_20D_U.xy")
    #Important Initial parameters for plottting
    VD_t = [.49, .48, .465, .36, .27, .23, .175, .165, .137]
    VD_pd = [.9, .87, .63, .37, .27, .235, .18, .17, .13]
    DDS = [2, 3, 4, 6, 8, 10, 15, 20, 25]
    U_inf = 6.82 #Free Stream velocity (m/s)
    D = .25 #Diameter of turbine (meters)
    #Labels to be reused 
    labelx = '$Diameters Down Stream y/D (D = 0.25m)$'
    labely = '$Velocity Defecit (U_\infty-U_c)/U_\infty$'
    labelk = '$Turbulent Kinetic Energy (m^2/s^2)$'
    labelw = '$Dissipation Rate (1/s)$'

    x, u, v, w = np.loadtxt(filepath, unpack=True)
    #u_nf = (U_inf-u)/U_inf
    #x_nf = x/D #full length of wind tunnel
    u_n = (U_inf-u[280:])/U_inf #after Turbine
    x_n = (x[280:]-35) / D #after Turbine

    plt.plot(x, u/U_inf, "r")
    plt.title("Velocity at Hub Height Along Wind Tunnel")
    plt.xlabel('Distance Along Tunnel (m)')
    plt.ylabel('$Normalized Velocity U/U_\infty$')
    plt.grid()
    plt.tight_layout()

    if name:
        plt.show()

    #plot of Velocity Defecit after Turbine
    plt.plot(DDS, VD_t, "-*b")
    plt.plot(DDS, VD_pd,'-om')
    plt.plot(x_n, u_n, "k")
    plt.title("Recovery of Velocity Deficit behind Actuator Disk")
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.grid()
    plt.tight_layout()

    if name:
        plt.show()

    #plot of TKE after Turbine
    x2, k, omega, p = np.loadtxt(filepath2, unpack=True)
    plt.plot(x_n, k[280:], "r")
    plt.title("Turbulent Kinetic Energy")
    plt.xlabel(labelx)
    plt.ylabel(labelk)
    plt.grid()
    plt.tight_layout()

    if name:
        plt.show()

    #plot of Profile behind Turbine
    h1, u1, v1, w1 = np.loadtxt(filepath3, unpack=True)
    h2, u2, v2, w2 = np.loadtxt(filepath4, unpack=True)
    plt.plot( u1/U_inf ,h1, "r")
    plt.plot( u2/U_inf ,h2, "g")
    plt.title("Velocity Profile behind Turbine")
    plt.xlabel('$Normalized Velocity U/U_\infty$')
    plt.ylabel("Height (m)")
    plt.grid()
    plt.tight_layout()

    if name:
        plt.show()



def plot_velocity():
    plt.figure()
    x = np.linspace(0, 100)
    y = np.sin(x)
    plt.plot(x, y)
    plt.xlabel("Xlabel")
    plt.show()
    
def something():
    """This is a function."""
    print("Something")

if __name__ == "__main__":
    load_text_file("hub_U.xy")
    #plot_velocity()
    #something()
