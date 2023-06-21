import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint, quad

def dydt(y, t, omg, d):
    # \dot\dot x = - omg^2 * x - d * \dot x 

    # system of 1. order ode
    # \dot x = y
    # \dot y = - \omg^2 x - d * y

    return [ y[1], - omg**2 * y[0] - d * y[1] ]

t_str = 0. # starting
t_end = 6. # end time

dt = (t_end - t_str)/1000. # time increment for output

# time array for output
tt = np.arange(t_str, t_end, dt)

omg = 6. # frequency omega
d = 1.   # damping constant     

# initial condition
#      x(0) = 0
# \dot x(0) = 1

y0 = [1, 0.]

y = odeint(dydt, y0, tt, args=(omg, d))

fig, ax = plt.subplots( 2 , 1 )

ax[0].plot(tt, y[:, 0])

ax[1].plot(y[:, 0], y[:, 1])

plt.savefig('ex10.pdf')
