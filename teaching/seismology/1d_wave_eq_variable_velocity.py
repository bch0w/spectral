""" 
See UAFGEOTEACH/GEOS604_seismo_soln/solutions/geos604funcs for the most up to 
date version of this

This was copied directly from: 
ttps://github.com/sachabinder/wave_equation_simulations/blob/main/1D_WAVE-EQ_variable-velocity.py

Citation: Sacha BINDER. Étude de l’observation et de la modélisation des ondes 
de surface en eau peu profonde. Physique-Informatique.TIPE session 2021.

I have modified it slightly to work behind the scenes with a Jupyter notebook 
for UAF GEOS604

__________________

This file was built to solve numerically a classical PDE, 1D wave equation. 

The numerical scheme is based on finite difference method. This program is also 
providing several boundary conditions. More particularly the Neumann, Dirichlet 
and Mur boundary conditions.

Copyright - © SACHA BINDER - 2021
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def anim_1D(x,y, pas_de_temps, pas_d_images, junction=None, save=False, 
            myxlim=(0, 4), myylim = (-4, 4)):
    """
    Function allowing to display an annimation based on calculation result with 
    a given time step. This function can be used to save the images sequence in 
    the current directory.
    
    The y parameter is a list containing several functions to display : y
     = [ [f_1(x)], ... , [f_n(x)] ].
    
    (x:np.ndarray (format 1D), y:np.ndarray (format 2D), pas_de_temps:float , 
    pas_d_images:int, save:bool , myxlim:tuple , myylim:tuple) -> plot (+ .mp4)
    """
    fig = plt.figure()
    ax = plt.axes(xlim= myxlim , ylim= myylim)
    line, = ax.plot([], [])
    ax.set_title("t = 0 s")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("$u$ [m]")
    if junction:
        plt.axvline(junction, ls="--", c="k")

    def init():
        line.set_data([],[])
        return line,
    
    # animation function.  This is called sequentially
    def animate(i):
        line.set_data(x, y[:,pas_d_images*i])
        ax.set_title("t = {:.2f} s".format(np.round(i*pas_d_images*pas_de_temps, 4)))
        return line,
        
    # call the animator. blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, 
                                   frames=y.shape[1]//pas_d_images, interval=10)

    return anim

def gaussian(x, b):
    """
    Single space variable fonction that 
    represent the wave form at t = 0 with starting location `b`
    """
    return np.exp(-(x-b)**2/0.01)

def plot_1d_wave_eq(length_of_string=1.5, duration_of_simulation=4,
                    junction=0.7, c_1=1.0, c_2=0.5, source_location=0,
                    source_location_2=0,
                    left_bound_cond="neumann", right_bound_cond="dirichlet"):
    """
    Solve the 1D Wave equation and plot
    """
    loop_exec = 1  # Processing loop execution flag

    assert(left_bound_cond in ["neumann", "dirichlet", "mur"])
    assert(right_bound_cond in ["neumann", "dirichlet", "mur"])

    #Spatial mesh - i indices
    L_x = length_of_string #Range of the domain according to x [m]
    dx = 0.01 #Infinitesimal distance
    N_x = int(L_x/dx) #Points number of the spatial mesh
    X = np.linspace(0,L_x,N_x+1) #Spatial array

    #Temporal mesh with CFL < 1 - j indices
    L_t = duration_of_simulation #Duration of simulation [s]
    dt = 0.01 * dx  # Infinitesimal time with CFL (Courant–Friedrichs–Lewy condition)
    N_t = int(L_t/dt) #Points number of the temporal mesh
    T = np.linspace(0,L_t,N_t+1) #Temporal array

    #Velocity array for calculation (finite elements)
    c = np.zeros(N_x+1, float)
    for i in range(0,N_x+1):
        if X[i] <= junction:
            c[i] = c_1
        else:
            c[i] = c_2

    ############## CALCULATION CONSTANTS ###############
    C2 = (dt/dx)**2

    CFL_1 = c_1*(dt/dx)
    CFL_2 = c_2*(dt/dx)
    # print(f"CFL_1={CFL_1}")
    # print(f"CFL_2={CFL_2}")

    ############## PROCESSING LOOP ###############

    if loop_exec:
        # $\forall i \in {0,...,N_x}$
        u_jm1 = np.zeros(N_x+1,float)   #Vector array u_i^{j-1}
        u_j = np.zeros(N_x+1,float)     #Vector array u_i^j
        u_jp1 = np.zeros(N_x+1,float)   #Vector array u_i^{j+1}
        
        q = np.zeros(N_x+1,float)
        q[0:N_x+1] = c[0:N_x+1]**2
        
        U = np.zeros((N_x+1,N_t+1),float) #Global solution
        
        #init cond - at t = 0
        u_j[0:N_x+1] = gaussian(X[0:N_x+1], b=source_location)

        U[:,0] = u_j.copy()
        
        
        #init cond - at t = 1
        #without boundary cond
        u_jp1[1:N_x] = u_j[1:N_x] + 0.5*C2*( 0.5*(q[1:N_x] + q[2:N_x+1])*(u_j[2:N_x+1] - u_j[1:N_x]) - 0.5*(q[0:N_x-1] + q[1:N_x])*(u_j[1:N_x] - u_j[0:N_x-1]))
        
        #left boundary conditions
        if left_bound_cond == "dirichlet":
            u_jp1[0] = 0
        elif left_bound_cond == "neumann":
            u_jp1[0] = u_j[0] + 0.5*C2*( 0.5*(q[0] + q[0+1])*(u_j[0+1] - u_j[0]) - 0.5*(q[0] + q[0+1])*(u_j[0] - u_j[0+1]))
        elif left_bound_cond == "mur":
            u_jp1[0] = u_j[1] + (CFL_1 -1)/(CFL_1 + 1)*( u_jp1[1] - u_j[0])
        
        #right boundary conditions
        if right_bound_cond == "dirichlet":
            u_jp1[N_x] = 0
        elif right_bound_cond == "neumann":
            u_jp1[N_x] =  u_j[N_x] + 0.5*C2*( 0.5*(q[N_x-1] + q[N_x])*(u_j[N_x-1] - u_j[N_x]) - 0.5*(q[N_x-1] + q[N_x])*(u_j[N_x] - u_j[i-1])) 
        elif right_bound_cond == "mur":
            u_jp1[N_x] = u_j[N_x-1] + (CFL_2 -1)/(CFL_2 + 1)*(u_jp1[N_x-1] - u_j[N_x])
        
        u_jm1 = u_j.copy()  #go to the next step
        u_j = u_jp1.copy()  #go to the next step
        U[:,1] = u_j.copy()
        
        #Process loop (on time mesh)
        for j in range(1, N_t):
            #calculation at step j+1
            #without boundary cond
            u_jp1[1:N_x] = -u_jm1[1:N_x] + 2*u_j[1:N_x] + C2*( 0.5*(q[1:N_x] + q[2:N_x+1])*(u_j[2:N_x+1] - u_j[1:N_x]) - 0.5*(q[0:N_x-1] + q[1:N_x])*(u_j[1:N_x] - u_j[0:N_x-1]))
            
            
            #left bound conditions
            if left_bound_cond == "dirichlet":
                #Dirichlet bound cond
                u_jp1[0] = 0
            elif left_bound_cond == "neumann":
                #Nuemann bound cond
                #i = 0
                u_jp1[0] = -u_jm1[0] + 2*u_j[0] + C2*( 0.5*(q[0] + q[0+1])*(u_j[0+1] - u_j[0]) - 0.5*(q[0] + q[0+1])*(u_j[0] - u_j[0+1]))       
            elif left_bound_cond == "mur":
                #Mur bound cond
                #i = 0
                u_jp1[0] = u_j[1] + (CFL_1 -1)/(CFL_1 + 1)*( u_jp1[1] - u_j[0])

            #right bound conditions
            if right_bound_cond == "dirichlet":
                #Dirichlet bound cond
                u_jp1[N_x] = 0               
            elif right_bound_cond == "neumann":
                #Nuemann bound cond
                #i = N_x
                u_jp1[N_x] = -u_jm1[N_x] + 2*u_j[N_x] + C2*( 0.5*(q[N_x-1] + q[N_x])*(u_j[N_x-1] - u_j[N_x]) - 0.5*(q[N_x-1] + q[N_x])*(u_j[N_x] - u_j[N_x-1]))            
            elif right_bound_cond == "mur":
                #Mur bound cond
                #i = N_x
                u_jp1[N_x] = u_j[N_x-1] + (CFL_2 -1)/(CFL_2 + 1)*(u_jp1[N_x-1] - u_j[N_x])
            
            u_jm1[:] = u_j.copy()   #go to the next step
            u_j[:] = u_jp1.copy()   #go to the next step
            U[:,j] = u_j.copy()

    # Plot
    anim1 = anim_1D(X, U, dt, 10, save=False, myxlim=(0, length_of_string), 
                    myylim=(-1.0,1.5), junction=junction)
    plt.show()


if __name__ == "__main__":
    plot_1d_wave_eq()
