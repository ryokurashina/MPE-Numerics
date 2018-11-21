""""" File contains functions that produces plots for the report """""

# Import packages and libraries
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

from initialConditions import *
from linAdSchemes import *
from sWSchemes import *
from analyse import *
from math import pi

def dispersion(N, k, dt, n_steps, H):
    """ Produces plots to be able to visibly see dispersion errors """
    # Set grid and parameters
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Initialise grid-size and total time
    dx = 1/(N-1)
    dt = 1e-3
    n_steps = 100
    T = n_steps*dt

    x_shift = x + 0.5*dx*np.ones_like(x)

    # Courant number
    c = np.sqrt(g*H)*dt/dx
    print("Courant Number (c): %f" %(c))

    # SW
    u1, h1 = UFB(u_init, h_init, n_steps, c, H)
    u2, h2 = SFB(u_init, h_init, n_steps, c, H, x, x_shift)
    u_exact, h_exact = trav_wave(x, T, k, H, 1)

    # Plot results
    plt.figure(1)
    plt.plot(x,u1, label="UFB")
    plt.plot(x,u2, label="SFB")
    plt.plot(x,u_exact, label="Exact Solution")
    plt.plot(x,u_init, '--', label="Initial Condition")
    plt.title('(a)')
    plt.xlabel('$x$')
    plt.ylabel('$u$')
    plt.legend()
    plt.savefig("plots/Dispersion_u.png")
    plt.close()

    plt.figure(2)
    plt.plot(x, h1, label="UFB")
    plt.plot(x, h2, label="SFB")
    plt.plot(x ,h_exact, label="Exact Solution")
    plt.plot(x ,h_init, '--', label="Initial Condition")
    plt.title('(b)')
    plt.xlabel('$x$')
    plt.ylabel('$h$')
    plt.legend()
    plt.savefig("plots/Dispersion_h.png")
    plt.close()

def L1_growth(N, k, dt, n_steps, H):
    """ Produces plots which tracks the L1-norm error growth in time """
    # Set grid and parameters
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Initialise grid-size and total time
    dx = 1/(N-1)
    T = n_steps*dt

    x_shift = x + 0.5*dx*np.ones_like(x)

    # Courant number
    c = np.sqrt(g*H)*dt/dx
    print("Courant Number (c): %f" %(c))

    # Array to store L2-errors
    U1 = np.zeros(n_steps)
    U2 = np.zeros(n_steps)
    U_exact = np.zeros(n_steps)
    H1 = np.zeros(n_steps)
    H2 = np.zeros(n_steps)
    H_exact = np.zeros(n_steps)

    # Initialise for loop
    u_old1 = u_init
    h_old1 = h_init
    u_old2 = u_init
    h_old2 = h_init

    # Perform a UFB and SFB step each iteration
    for i in range(n_steps):
        t = (i+1)*dt
        # SW
        u_new1, h_new1 = UFB(u_old1, h_old1, 1, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, 1, c, H , x, x_shift)
        u_exact, h_exact = trav_wave(x, t, k, H, 1)

        # Store the L1-norms
        U1[i] = total_abs(u_new1)
        U2[i] = total_abs(u_new2)
        U_exact[i] = total_abs(u_exact)

        H1[i] = total_abs(h_new1)
        H2[i] = total_abs(h_new2)
        H_exact[i] = total_abs(h_exact)

        u_old1 = u_new1.copy()
        h_old1= h_new1.copy()

        u_old2 = u_new1.copy()
        h_old2= h_new1.copy()

    # Plot results
    plt.figure(1)
    plt.plot(range(1,n_steps+1), U1, label="UFB")
    plt.plot(range(1,n_steps+1), U2, label="SFB")
    plt.plot(range(1,n_steps+1), U_exact, "--", label="Exact Solution")
    plt.title('(a)')
    plt.xlabel('Time-step')
    plt.ylabel('$L^1(u)$')
    plt.legend()
    plt.savefig("plots/L1_growth_u.png")
    plt.close()

    # Plot results
    plt.figure(2)
    plt.plot(range(1,n_steps+1), H1, label="UFB")
    plt.plot(range(1,n_steps+1), H2, label="SFB")
    plt.plot(range(1,n_steps+1), H_exact, "--", label="Exact Solution")
    plt.title('(b)')
    plt.xlabel('Time-step')
    plt.ylabel('$L^1(h)$')
    plt.legend()
    plt.savefig("plots/L1_growth_h.png")
    plt.close()

def L2_growth(N, k, dt, n_steps, H):
    """
    Produces plots which tracks the L2-norm error growth in time
    Not used in final report
    """
    # Set grid and parameters
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81
    H = 1

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Initialise grid-size and total time
    dx = 1/(N-1)
    T = n_steps*dt

    x_shift = x + 0.5*dx*np.ones_like(x)

    # Courant number
    c = np.sqrt(g*H)*dt/dx
    print("Courant Number (c): %f" %(c))

    # Array to store L2-errors
    e1_u = np.zeros(n_steps)
    e2_u = np.zeros(n_steps)
    e1_h = np.zeros(n_steps)
    e2_h = np.zeros(n_steps)

    u_old1 = u_init
    h_old1 = h_init
    u_old2 = u_init
    h_old2 = h_init

    # Perform one UFB and SFB step per iteration
    for i in range(n_steps):
        t = (i+1)*dt
        # SW
        u_new1, h_new1 = UFB(u_old1, h_old1, 1, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, 1, c, H , x, x_shift)
        u_exact, h_exact = trav_wave(x, t, k, H, 1)

        # Store L2-norm values
        e1_u[i] = L2_error(u_new1, u_exact)
        e2_u[i] = L2_error(u_new2, u_exact)

        e1_h[i] = L2_error(h_new1, h_exact)
        e2_h[i] = L2_error(h_new2, h_exact)

        u_old1 = u_new1.copy()
        h_old1= h_new1.copy()

        u_old2 = u_new1.copy()
        h_old2= h_new1.copy()

    # Plot results
    plt.figure(1)
    plt.plot(range(1,n_steps+1), e1_u, label="UFB")
    plt.plot(range(1,n_steps+1), e2_u, label="SFB")
    plt.title('(a)')
    plt.xlabel('Time-step')
    plt.ylabel('$L^2(u-u_{exact})$')
    plt.legend()
    plt.close()

    # Plot results
    plt.figure(2)
    plt.plot(range(1,n_steps+1), e1_h, label="UFB")
    plt.plot(range(1,n_steps+1), e2_h, label="SFB")
    plt.title('(b)')
    plt.xlabel('Time-step')
    plt.ylabel('$L^2(h-h_{exact})$')
    plt.legend()
    plt.close()

def L2_convergence_x(k, dt, n_steps, H):
    """ Shows the L2-norm convergence of UFB and SFB schemes in space """
    # Set resolution and parameters
    N_vec = [8, 16, 32, 64, 128, 256, 512, 1024]
    M = len(N_vec)
    g = 9.81

    # Initialise timestep, grid-size and advection velocity
    T = n_steps*dt

    # Array to store L2-errors
    e1_u = np.zeros(M)
    e2_u = np.zeros(M)
    e1_h = np.zeros(M)
    e2_h = np.zeros(M)

    for i in range(M):
        x_ = np.linspace(0, 1, N_vec[i]+1)
        x = x_[0:N_vec[i]]
        dx = 1/(N_vec[i]-1)
        x_shift = x + 0.5*dx*np.ones_like(x)

        # Initial condition
        k = 4
        u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
        h_init = np.cos(2*pi*k*x)

        u_old1 = u_init
        h_old1 = h_init
        u_old2 = u_init
        h_old2 = h_init

        # Courant number
        c = np.sqrt(g*H)*dt/dx
        print("Courant Number (c): %f" %(c))

        # SW
        u_new1, h_new1 = UFB(u_old1, h_old1, n_steps, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, n_steps, c, H , x, x_shift)
        u_exact, h_exact = trav_wave(x, T, k, H, 1)

        e1_u[i] = L2_error(u_new1, u_exact)
        e2_u[i] = L2_error(u_new2, u_exact)

        e1_h[i] = L2_error(h_new1, h_exact)
        e2_h[i] = L2_error(h_new2, h_exact)

    # Plot results
    plt.figure(1)
    plt.plot(np.log(N_vec), np.log(e1_u), "^-", label="UFB")
    plt.plot(np.log(N_vec), np.log(e2_u), "^-", label="SFB")
    plt.plot(np.log(N_vec), -2*np.log(N_vec)+2, label="Line Gradient -2")
    plt.title('(a)')
    plt.xlabel('$\log(N)$')
    plt.ylabel('$L^2(u-u_{exact})$')
    plt.legend()
    plt.savefig("plots/spatial_conv_u")
    plt.close()

    # Plot results
    plt.figure(2)
    plt.plot(np.log(N_vec), np.log(e1_h), "^-", label="UFB")
    plt.plot(np.log(N_vec), np.log(e2_h), "^-", label="SFB")
    plt.plot(np.log(N_vec), -2*np.log(N_vec), label="Line Gradient -2")
    plt.title('(b)')
    plt.xlabel('$\log(N)$')
    plt.ylabel('$L^2(h-h_{exact})$')
    plt.legend()
    plt.savefig("plots/spatial_conv_h")
    plt.close()

def L2_convergence_t(N, k, n_steps, H):
    """
    Shows the L2-norm convergence of UFB and SFB schemes in time
    Not used in final report
    """
    # Set time-steps and parameters
    dt = 1e-4
    dt_vec = np.array([dt, dt/2, dt/4, dt/8, dt/16])
    M = len(dt_vec)
    g = 9.81

    # Array to store L2-errors
    e1_u = np.zeros(M)
    e2_u = np.zeros(M)
    e1_h = np.zeros(M)
    e2_h = np.zeros(M)

    # Set grid
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    dx = 1/(N-1)
    x_shift = x + 0.5*dx*np.ones_like(x)

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Perform UFB and SFB for each grid
    for i in range(M):
        # We want to get to the same point in time
        n_steps = int(2e-2/dt_vec[i])
        T = n_steps*dt_vec[i]

        u_old1 = u_init
        h_old1 = h_init
        u_old2 = u_init
        h_old2 = h_init

        # Courant number
        c = np.sqrt(g*H)*dt_vec[i]/dx
        print("Courant Number (c): %f" %(c))

        # SW
        u_new1, h_new1 = UFB(u_old1, h_old1, n_steps, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, n_steps, c, H , x, x_shift)
        u_exact, h_exact = trav_wave(x, T, k, H, 1)

        # Store L2-norm errors for each grid
        e1_u[i] = L2_error(u_new1, u_exact)
        e2_u[i] = L2_error(u_new2, u_exact)

        e1_h[i] = L2_error(h_new1, h_exact)
        e2_h[i] = L2_error(h_new2, h_exact)

    # Plot results
    plt.figure(1)
    plt.plot(np.log(dt_vec), np.log(e1_u), '^-', label="UFB")
    plt.plot(np.log(dt_vec), np.log(e2_u), '^-', label="SFB")
    plt.title('(a)')
    plt.xlabel('$\log(dt)$')
    plt.ylabel('$log(L^2(u-u_{exact}))$')
    plt.legend()
    plt.show()

    # Plot results
    plt.figure(2)
    plt.plot(np.log(dt_vec), np.log(e1_h), '^-', label="UFB")
    plt.plot(np.log(dt_vec), np.log(e2_h), '^-', label="SFB")
    plt.title('(b)')
    plt.xlabel('\log(dt)')
    plt.ylabel('$log(L^2(h-h_{exact}))$')
    plt.legend()
    plt.show()

def calc_moments(N, k, dt, n_steps, H, p):
    """ Calculates the centered p-th moment for each time-step and produces plots"""
    # Set grid and parameters
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Initialise grid-size and total time
    dx = 1/(N-1)
    T = n_steps*dt

    x_shift = x + 0.5*dx*np.ones_like(x)

    # Courant number
    c = np.sqrt(g*H)*dt/dx
    print("Courant Number (c): %f" %(c))

    # Array to store moments at each timestep
    m1_u = np.zeros(n_steps)
    m2_u = np.zeros(n_steps)
    m_u_exact = np.zeros(n_steps)
    m1_h = np.zeros(n_steps)
    m2_h = np.zeros(n_steps)
    m_h_exact = np.zeros(n_steps)

    u_old1 = u_init
    h_old1 = h_init
    u_old2 = u_init
    h_old2 = h_init

    # Calculate the centered p-th moment at each time-step for UFB and SFB
    for i in range(n_steps):
        t = (i+1)*dt
        # SW
        u_new1, h_new1 = UFB(u_old1, h_old1, 1, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, 1, c, H , x, x_shift)
        u_exact, h_exact = trav_wave(x, t, k, H, 1)

        m1_u[i] = p_moment(u_new1, p)
        m2_u[i] = p_moment(u_new2, p)
        m_u_exact[i] = p_moment(u_exact, p)

        m1_h[i] = p_moment(h_new1, p)
        m2_h[i] = p_moment(h_new2, p)
        m_h_exact[i] = p_moment(h_exact, p)

        u_old1 = u_new1.copy()
        h_old1= h_new1.copy()

        u_old2 = u_new1.copy()
        h_old2= h_new1.copy()

    # Plot results
    n_steps_vec = range(1,n_steps+1)

    plt.figure(1)
    plt.plot(n_steps_vec, m1_u, label="UFB")
    plt.plot(n_steps_vec, m2_u, label="SFB")
    plt.plot(n_steps_vec, m_u_exact, "--", label="Exact")
    plt.title('(a)')
    plt.xlabel('Time-step')
    plt.ylabel('p-th Moment (p = %s)' %(p))
    plt.legend()
    plt.savefig("plots/Moment_u.png")
    plt.close()

    plt.figure(2)
    plt.plot(n_steps_vec, m1_h, label="UFB")
    plt.plot(n_steps_vec, m2_h, label="SFB")
    plt.plot(n_steps_vec, m_h_exact, "--", label="Exact")
    plt.title('(b)')
    plt.xlabel('Time-step')
    plt.ylabel('p-th Moment (p = %s)' %(p))
    plt.legend()
    plt.savefig("plots/Moment_h.png")
    plt.close()
