"""""
File contains functions that produces plots for the report. Includes plots for:
* Dispersion properties of UFB and SFB
* L1-norm growth of UFB and SFB
* L2-norm growth of UFB and SFB
* L2-norm convergence of UFB and SFB with decreasing spatial-steps
* L2-norm convergence of UFB and SFB with decreasing time-steps
* Calculating the centered p-th moment of UFB and SFB schemes
"""""

# Import packages and libraries
import numpy as np
import matplotlib.pyplot as plt

from initialConditions import *
from linAdSchemes import *
from sWSchemes import *
from analyse import *
from math import pi

def dispersion(N, k, dt, n_steps, H):
    # Set resolution and parameters
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81
    H = 1

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Initialise timestep, grid-size and advection velocity
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
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend()
    plt.show()

    plt.figure(2)
    plt.plot(x, h1, label="UFB")
    plt.plot(x, h2, label="SFB")
    plt.plot(x ,h_exact, label="Exact Solution")
    plt.plot(x ,h_init, '--', label="Initial Condition")
    plt.xlabel('x')
    plt.ylabel('h')
    plt.legend()
    plt.show()

def L1_growth(N, k, dt, n_steps, H):
    # Set resolution and parameters
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Initialise timestep, grid-size and advection velocity
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

    u_old1 = u_init
    h_old1 = h_init
    u_old2 = u_init
    h_old2 = h_init

    for i in range(n_steps):
        t = (i+1)*dt
        # SW
        u_new1, h_new1 = UFB(u_old1, h_old1, 1, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, 1, c, H , x, x_shift)
        u_exact, h_exact = trav_wave(x, t, k, H, 1)

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
    plt.title('L1-norms for u')
    plt.xlabel('Timestep')
    plt.ylabel('L1-norm')
    plt.legend()
    plt.show()

    # Plot results
    plt.figure(2)
    plt.plot(range(1,n_steps+1), H1, label="UFB")
    plt.plot(range(1,n_steps+1), H2, label="SFB")
    plt.plot(range(1,n_steps+1), H_exact, "--", label="Exact Solution")
    plt.title('L1-norms for h')
    plt.xlabel('Timestep')
    plt.ylabel('L1-norm')
    plt.legend()
    plt.show()

def L2_growth(N, k, dt, n_steps, H):
    # Set resolution and parameters
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81
    H = 1

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Initialise timestep, grid-size and advection velocity
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

    for i in range(n_steps):
        t = (i+1)*dt
        # SW
        u_new1, h_new1 = UFB(u_old1, h_old1, 1, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, 1, c, H , x, x_shift)
        u_exact, h_exact = trav_wave(x, t, k, H, 1)

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
    plt.title('L2-error norms for u')
    plt.xlabel('Timestep')
    plt.ylabel('L2_error')
    plt.legend()
    plt.show()

    # Plot results
    plt.figure(2)
    plt.plot(range(1,n_steps+1), e1_h, label="UFB")
    plt.plot(range(1,n_steps+1), e2_h, label="SFB")
    plt.title('L2-error norms for h')
    plt.xlabel('Timestep')
    plt.ylabel('L2_error')
    plt.legend()
    plt.show()

def L2_convergence_x(k, dt, n_steps, H):
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
    plt.plot(np.log(N_vec), np.log(e1_u), label="UFB")
    plt.plot(np.log(N_vec), np.log(e2_u), label="SFB")
    plt.plot(np.log(N_vec), -2*np.log(N_vec)+2, label="Line Gradient -2")
    plt.title('L2-error norms for h')
    plt.xlabel('log(N)')
    plt.ylabel('log(L2-error)')
    plt.legend()
    plt.show()

    # Plot results
    plt.figure(2)
    plt.plot(np.log(N_vec), np.log(e1_h), label="UFB")
    plt.plot(np.log(N_vec), np.log(e2_h), label="SFB")
    plt.plot(np.log(N_vec), -2*np.log(N_vec), label="Line Gradient -2")
    plt.title('L2-error norms for h')
    plt.xlabel('log(N)')
    plt.ylabel('log(L2-error)')
    plt.legend()
    plt.show()

def L2_convergence_t(N, k, n_steps, H):
    # Set resolution and parameters
    dt = 1e-4
    dt_vec = np.array([dt, dt/2, dt/4, dt/8, dt/16])
    M = len(dt_vec)
    g = 9.81
    H = 1

    # Initialise timestep, grid-size and advection velocity
    # Array to store L2-errors
    e1_u = np.zeros(M)
    e2_u = np.zeros(M)
    e1_h = np.zeros(M)
    e2_h = np.zeros(M)

    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    dx = 1/(N-1)
    x_shift = x + 0.5*dx*np.ones_like(x)

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

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

        e1_u[i] = L2_error(u_new1, u_exact)
        e2_u[i] = L2_error(u_new2, u_exact)

        e1_h[i] = L2_error(h_new1, h_exact)
        e2_h[i] = L2_error(h_new2, h_exact)

    # Plot results
    plt.figure(1)
    plt.plot(np.log(dt_vec), np.log(e1_u), '^', label="UFB")
    plt.plot(np.log(dt_vec), np.log(e2_u), '^', label="SFB")
    # plt.plot(np.log(dt_vec), -2*np.log(dt_vec)+2, label="Line Gradient -2")
    plt.title('L2-error norms for u')
    plt.xlabel('log(dt)')
    plt.ylabel('log(L2-error)')
    plt.legend()
    plt.show()

    # Plot results
    plt.figure(2)
    plt.plot(np.log(dt_vec), np.log(e1_h), label="UFB")
    plt.plot(np.log(dt_vec), np.log(e2_h), label="SFB")
    # plt.plot(np.log(dt_vec), -2*np.log(dt_vec), label="Line Gradient -2")
    plt.title('L2-error norms for h')
    plt.xlabel('log(dt)')
    plt.ylabel('log(L2-error)')
    plt.legend()
    plt.show()

def calc_moments(N, k, dt, n_steps, H, p):
    # Set resolution and parameters
    # Calculate which moment
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81

    # Initial condition
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Initialise timestep, grid-size and advection velocity
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
    plt.plot(n_steps_vec, m_u_exact, label="Exact")
    plt.title('p-th moment of u')
    plt.xlabel('timestep')
    plt.ylabel('moment')
    plt.legend()
    plt.show()

    plt.figure(2)
    plt.plot(n_steps_vec, m1_h, label="UFB")
    plt.plot(n_steps_vec, m2_h, label="SFB")
    plt.plot(n_steps_vec, m_h_exact, label="Exact")
    plt.title('p-th moment of h')
    plt.xlabel('timestep')
    plt.ylabel('moment')
    plt.legend()
    plt.show()
