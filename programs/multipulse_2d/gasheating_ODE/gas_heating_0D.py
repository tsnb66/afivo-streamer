import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


data =  np.loadtxt("powerDensity_power_density_pointVsTime.txt")


def f(y, t, eta_f, eta_s, tau_vt):
    E, Evt = y
    jE = np.interp(t, data[:,0], data[:,1])
    return [eta_f*jE + Evt, eta_s*jE - Evt/tau_vt]




dt = 1e-10
tstart, tend = 0, 5e-6
Nt = int((tend-tstart)/dt)
t = np.linspace(tstart, tend, Nt)
eta_f, eta_s, tau_vt = 0.3, 0.7, 20e-6
E0 = 1e5/(1.4-1.0)


solution = odeint(f, [E0, 0], t, args=(eta_f, eta_s, tau_vt)) 


plt.plot(t, solution[:,0], 'b-', label="Translation")
plt.plot(t, solution[:,1], 'g-', label="Vibration")
plt.legend()
plt.figure(2)
kb = 1.38e-23
N = 2.4e25
Tg = (solution[:,0]*(1.4-1))/(kb*N)
P = solution[:,0]*0.4
plt.plot(t, Tg)
plt.show()





