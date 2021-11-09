# Maneuver 1

"""
Escape Maneuver: Exit Earth's Orbit
V_initial = v_0 = Orbital Velocity
V_final >= v_f = Escape Velocity
a_initial = F/m
Thrust = Constant = T
Specific Impulse = Constant = I_sp
Mass of Earth = M_earth = 5.972 × 10^24 kg

Orbital Period = T = 90 mins
Exhaust Velocity = v_e = g_0 I_sp
Rate of Change of Mass = p = dm/dt = T/v_e
v_0 - v_f = v_e * ln(m_0/m_f)
Taking Excape Velocity = v_f = 11.2 km/s
Acceleration = a = - p v_e / m(t)
"""

from astropy import units as u
import numpy as np

#----------------------------------------------------
g_0 = 9.80665 * u.m / (u.second * u.second)
G = 6.67408e-11 * u.m**3 * u.kg**(-1) * u.s**(-2)
M_Earth = 5.972e24 * u.kg
R_Earth = (6371 * u.km).to(u.m)
#----------------------------------------------------

#----------------------------------------------------
TP_orbit = 29499.67820209 * u.s # (90 * u.min).to(u.s)
m_0 = 5000 * u.kg
I_sp = 4000 * u.second
T = 400e-3 * u.kg * u.m  / (u.second * u.second) # newtons
v_f = 11.2e3 * u.m / u.second
#---------------------------------------------------

R_orbit = (TP_orbit*TP_orbit*G*M_Earth / (4*np.pi*np.pi))**(1/3)

v_0 = np.sqrt(G*M_Earth/R_orbit)

v_e = g_0 * I_sp

p = T / v_e

m_f = m_0 * np.exp(-1 * (v_f - v_0) / v_e)

Time_Elapsed = (m_0 - m_f)/p

a_0 = -1 * p * v_e/ m_0
a_f = -1 * p * v_e/ m_f

print(m_0, m_f)
print(a_0, a_f)
print(Time_Elapsed)