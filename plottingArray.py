#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import math
import scipy.io


def plot_vel_def():
    #This Function loads and plots specified result data
    name = "hub_U.xy"
    print("Loading text file", name)
    #Read in data and sort it to get the Converged data file
    times = os.listdir("postProcessing/sets")
    times.sort()
    latest = times[-1]
    #Select the path to the specific file
    filepath = os.path.join("postProcessing", "sets", latest, 
                            name)
    #Important Initial parameters for plottting
    #read in defecit data
    FSfile1 = "/media/tom/b7e5b88c-2a56-43b2-8ed7-24c8be25d90f/DATAPlotting/FS_PD"
    FSpd = scipy.io.loadmat(FSfile1)
    FS_PD = FSpd['hub_def_FSPD']
    FS_PD = FS_PD[0]
    Diam_down = FSpd['Diam_down']
    Diam_down = Diam_down[0]
    FSfile2 = "/media/tom/b7e5b88c-2a56-43b2-8ed7-24c8be25d90f/DATAPlotting/FS_Turb"
    FSt = scipy.io.loadmat(FSfile2)
    FS_T = FSt['hub_def_FST']
    FS_T = FS_T[0]
    BLfile1 = "/media/tom/b7e5b88c-2a56-43b2-8ed7-24c8be25d90f/DATAPlotting/BL_PD"
    BLpd = scipy.io.loadmat(BLfile1)
    BL_PD = BLpd['hub_def_BLPD']
    BL_PD = BL_PD[0]
    BLfile2 = "/media/tom/b7e5b88c-2a56-43b2-8ed7-24c8be25d90f/DATAPlotting/BL_Turb"
    BLt = scipy.io.loadmat(BLfile2)
    BL_T = BLt['hub_def_BLT']
    BL_T = BL_T[0]
    #h_pd = pd['loc']
    #h_pd = h_pd[0]
                            
    VD_t = [.49, .48, .465, .36, .27, .23, .175, .165, .137]
    VD_pd = [.9, .87, .63, .37, .27, .235, .18, .17, .13]
    DDS = [2, 3, 4, 6, 8, 10, 15, 20, 25]
    D = .25 #Diameter of turbine (meters)
    turb_x = 39 #x-coordinate of the turbine
    #Labels to be reused 
    labelx = '$Diameters Down Stream y/D (D = 0.25m)$'
    labely = '$Velocity Defecit (U_\infty-U_c)/U_\infty$'
    #labelk = '$Turbulent Kinetic Energy (m^2/s^2)$'
    #labelw = '$Dissipation Rate (1/s)$'

    #load data
    x, u, v, w = np.loadtxt(filepath, unpack=True)
    
    #Find start of the turbine
    for i in range(1, len(u)):
        if x[i] == turb_x:
            #print "Equaled %d " % (turb_x)
            break
        elif x[i] > turb_x:
            #print "We're  greater than %d" % (turb_x)
            break
    U_inf = u[i-1] #Incoming velocity (m/s)
    start = i 
    u_n = (U_inf-u[start:])/U_inf #after Turbine
    x_n = (x[start:]-turb_x) / D #after Turbine

    #plot of Velocity Defecit after Turbine
    plt.figure()
    Turb, = plt.plot(DDS, VD_t, "b")
    #TurbBL, = plt.plot(Diam_down, BL_T, "->b")
    #TurbFS, = plt.plot(Diam_down, FS_T, "-*b")
    PD, = plt.plot(DDS, VD_pd,'r')
    #PDBL, = plt.plot(Diam_down, BL_PD,'->m')
    #PDFS, = plt.plot(Diam_down, FS_PD,'-*m')
    Hub, = plt.plot(x_n, u_n, "-+k")
    plt.title("Array: Recovery of Velocity Deficit 3rd Turbine at Hub Height", size=20)
    plt.xlabel(labelx, size=20)
    plt.ylabel(labely, size=20)
    plt.grid()
    plt.tight_layout()
    plt.legend( (Turb, PD, Hub), ('Model Turbine','Porous Disk','Numerical Model'),'best') #plt.legend( (Turb,TurbBL, TurbFS, PD,PDBL,PDFS, Hub), ('Turbine','Turbine BL','Turbine FS', 'Porous Disk','Porous Disk BL','Porous Disk FS','Simulation Results'),'best')

    if name:
        plt.show()

def plot_tunnel_plots():
    name = "hub_U.xy"
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
    #Important Initial parameters for plottting
    turb_x = 35 #x-coordinate of the turbine
    D = .25 #Diameter of turbine (meters)
    #Labels to be reused 
    labelx = r'Diameters Downstream $y/D$ ($D = 0.25$ m)'
    #labely = r'Velocity Deficit $(U_\infty-U_c)/U_\infty$'
    labelk = r'Turbulent Kinetic Energy $(\mathrm{m}^2/\mathrm{s}^2)$'
    #labelw = r'Dissipation Rate $(1/\mathrm{s})$'

    #Plot Velocity at Hub height in wind tunnel 
    x, u, v, w = np.loadtxt(filepath, unpack=True)
    U_inf = u[1] #Incoming velocity (m/s)
    for i in range(1, len(u)):
        if x[i] == turb_x:
            #print "Equaled %d " % (turb_x)
            break
        elif x[i] > turb_x:
            #print "We're  greater than %d" % (turb_x)
            break
    start = i 
    #u_n = (U_inf-u[start:])/U_inf #after Turbine
    x_n = (x[start:]-35) / D #after Turbine
    
    plt.figure()
    Hub, = plt.plot(x, u/U_inf, "r", label='Velocity at Hub Height')
    plt.title("Velocity at Hub Height Along Wind Tunnel")
    plt.xlabel('Distance Along Tunnel (m)')
    plt.ylabel(r'Normalized Velocity $U/U_\infty$')
    plt.grid()
    plt.tight_layout()
    #plt.legend( handles = [Hub])
    if name:
        plt.show()

    #plot of TKE after Turbine
    plt.figure()
    x_e, epsilon, k, p = np.loadtxt(filepath2, unpack=True)
    plt.plot(x_n, k[start:], "r")
    plt.title("Turbulent Kinetic Energy")
    plt.xlabel(labelx)
    plt.ylabel(labelk)
    plt.grid()
    plt.tight_layout()
    #plt.legend( (kHub), ('TKE at Hub Height'),'best')

    if name:
        plt.show()

def plot_profiles(Diam):
    #This Function loads and plots specified result data
    print("Loading text files")
    #Read in data and sort it to get the Converged data file
    times = os.listdir("postProcessing/sets")
    times.sort()
    latest = times[-1]
    #Select the path to the specific file
    sn = ("data_","_U.xy")
    name = Diam.join(sn)
    filepath1 = os.path.join("postProcessing", "sets", latest, 
                            name)
    st = ("/media/tom/b7e5b88c-2a56-43b2-8ed7-24c8be25d90f/DATAPlotting/Turbine Arrays/velocity_fullarray_8D_5x5_3rd_x",".mat")
    fileturb = Diam.join(st)
    #sp = ("/media/tom/b7e5b88c-2a56-43b2-8ed7-24c8be25d90f/DATAPlotting/Porous Disk BL/normvel_Disk_",".mat")
    #filepd = Diam.join(sp)
    fileinc = ("/media/tom/b7e5b88c-2a56-43b2-8ed7-24c8be25d90f/DATAPlotting/Porous Disk BL/normvel_incoming.mat")
    
   # turb_x = 35 #x-coordinate of the turbine
   #D = .25 #Diameter of turbine (meters)
    #Labels to be reused 
    lb = ("Array: Velocity Profile "," behind Turbine (Hub Height)")
    label_plot = Diam.join(lb)
    U_inf = 6.82
    #read in case resutls
    h, u, v, w = np.loadtxt(filepath1, unpack=True)
    #read in profile data
    t = scipy.io.loadmat(fileturb)
    v_t = t['meanvel']
    v_t = v_t[0]
    h_t = t['loc']
    h_t = h_t[0]
    #pd = scipy.io.loadmat(filepd)
    #v_pd = pd['NORMvel']
    #v_pd = v_pd[0]
    #h_pd = pd['loc']
    #h_pd = h_pd[0]
    inc = scipy.io.loadmat(fileinc)
    v_in = inc['NORMvel']
    v_in = v_in[0]
    h_in = inc['loc']
    h_in = h_in[0]
    lastt =len(v_t)
    
    u_plot =u[6:]/U_inf
    h_plot = h[6:]-.1875
    #print h_plot[1:]
    #print u_plot[1:]
    
    v_adj = v_t/v_t[lastt-1]
    h_adj = h_t-.1875
    #plot of Velocity Profile behind Turbine
    plt.figure()
    Sim, = plt.plot( u_plot ,h_plot, "-*k")
    Turb, = plt.plot(v_adj, h_adj, "r")
    #PD, = plt.plot(v_pd, h_pd, 'r')
    #plt.plot(v_in, h_in, 'r')
    plt.title(label_plot,size=20)
    plt.xlabel(r'Normalized Velocity $U/U_\infty$', size=20)
    plt.ylabel("$h/D (From Hub Height)$", size=20)
    plt.grid()
    plt.tight_layout()
    plt.legend( (Turb, Sim), ('3rd Model Turbine','Numerical Model'),'best')
    plt.show()
    
    low = 1000
    v_m= [0, 0 ,0, 0 ,0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    h_m= [0, 0 ,0, 0 ,0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #v_m= [0, 0 ,0, 0 ,0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]#uncomment to plot 4D
    #h_m= [0, 0 ,0, 0 ,0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]#uncomment to plot 4D
    #print u_plot
    #print h_plot
    #print v_m
    #print h_m
    #print h_adj[1]
    #Match data
    for i in range(1, len(h_plot)):
        low = 1000
        #print h_plot
        for j in range(1, len(h_adj)):
            
            if (h_plot[i]-h_adj[j])*(h_plot[i]-h_adj[j]) < low:
                
                
                low = (h_plot[i]-h_adj[j])*(h_plot[i]-h_adj[j])
                v_m[i]= v_adj[j]
                h_m[i]= h_adj[j]
                #print u_plot[i]
                #print v_m[i]
    #print h_plot
    #print v_n[46]
    #print v_m
    #print h_m
    #print len(v_m)
    #print len(h_m)     
    diff = 100*abs(v_m[1:]-u_plot[1:])/v_m[1:]
    plt.figure()
    #plt.plot( u_plot ,h_plot, "-*k")
    plt.plot( diff ,h_m[1:], "-*r")
    plt.title('Array: 3rd Turbine Profile Error (2D, Hub Height)',size=20)
    plt.xlabel('Velocity Percent Error (%)', size=20)
    plt.ylabel("$h/D (From Hub Height)$", size=20)
    plt.grid()
    #plt.axes([0, 10, 0, 1])
    plt.xlim(0, 10)
    plt.ylim(0, .9)
    plt.tight_layout()
    plt.show()
    #plot of Incoming Velocity Profile behind Turbine
    #plt.figure()
    #plt.plot( v_in, h_in, "b")
    #plt.title('BL: Incoming Velocity',size=20)
    #plt.xlabel(r'Normalized Velocity $U/U_\infty$', size=20)
    #plt.ylabel("$h/D  (D = 0.25m)$", size=20)
    #plt.grid()
    #plt.tight_layout()
    #plt.show()

if __name__ == "__main__":
    
  
    plot_vel_def()
    
    #plot_tunnel_plots()
    plot_profiles("2D")
    #plot_profiles("4D")
