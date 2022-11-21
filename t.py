import numpy as np
# Physical parameters
alpha = 0.1                     # Heat transfer coefficient
lx = 1.                         # Size of computational domain

# Grid parameters
nx = 21                         # number of grid points 
dx = lx / (nx-1)                # grid spacing
x = np.linspace(0., lx, nx)     # coordinates of grid points

# Time parameters
ti = 0.                         # initial time
tf = 5.                         # final time
fourier = 0.49                  # Fourier number
dt = fourier*dx**2/alpha        # time step
nt = int((tf-ti) / dt)          # number of time steps

# Initial condition
T0 = np.sin(2*np.pi*x)          # initial condition
source = 2*np.sin(np.pi*x)      # heat source term


def exact_solution(x,t, alpha):
    """Returns the exact solution of the 1D
    heat equation with heat source term sin(np.pi*x)
    and initial condition sin(2*np.pi*x)
    
    Parameters
    ----------
    x : array of floats
        grid points coordinates
    t: float
        time
    
    Returns
    -------
    f : array of floats
        exact solution
    """
    # Note the 'Pythonic' way to break the long line. You could
    # split a long line using a backlash (\) but the conventional
    # way is to embrace your code in parenthesis.
    #
    # For more info we refer you to PEP8:
    # https://www.python.org/dev/peps/pep-0008/#id19 
    f = (np.exp(-4*np.pi**2*alpha*t) * np.sin(2*np.pi*x)
      + 2.0*(1-np.exp(-np.pi**2*alpha*t)) * np.sin(np.pi*x) 
      / (np.pi**2*alpha))
    
    return f

print(exact_solution(x, 1, 1))
