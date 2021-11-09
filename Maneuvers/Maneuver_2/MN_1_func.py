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

def get_cost(m_0, TP_orbit, v_f):

    #----------------------------------------------------
    g_0 = 9.80665 * u.m / (u.second * u.second)
    G = 6.67408e-11 * u.m**3 * u.kg**(-1) * u.s**(-2)
    M_Earth = 5.972e24 * u.kg
    R_Earth = (6371 * u.km).to(u.m)
    #----------------------------------------------------

    #----------------------------------------------------
    # TP_orbit = 29499.67820209 * u.s # (90 * u.min).to(u.s)
    # m_0 = 5000 * u.kg
    TP_orbit = TP_orbit.to(u.s)
    m_0 = m_0.to(u.kg)
    v_f = v_f.to(u.m/u.s)
    
    I_sp = 4000 * u.second
    T = 400e-3 * u.kg * u.m  / (u.second * u.second) # newtons
    # v_f = 11.2e3 * u.m / u.second
    #---------------------------------------------------

    R_orbit = (TP_orbit*TP_orbit*G*M_Earth / (4*np.pi*np.pi))**(1/3)

    v_0 = np.sqrt(G*M_Earth/R_orbit)

    v_e = g_0 * I_sp

    p = -1 * T / v_e

    m_f = m_0 * np.exp(-1 * (v_f - v_0) / v_e)

    Time_Elapsed = (m_f - m_0)/p

    a_0 = -1 * p * v_e/ m_0
    a_f = -1 * p * v_e/ m_f

    print(m_0, m_f)
    print(a_0, a_f)
    print(v_0, v_f)
    print(Time_Elapsed)

# get_cost(5000*u.kg, 90*u.min, 19.5*u.km/u.s)

def get_cost_by_alt(m_0, TP_orbit, to_alt):

    #----------------------------------------------------
    g_0 = 9.80665 * u.m / (u.second * u.second)
    G = 6.67408e-11 * u.m**3 * u.kg**(-1) * u.s**(-2)
    M_Earth = 5.972e24 * u.kg
    R_Earth = (6371 * u.km).to(u.m)
    #----------------------------------------------------

    #----------------------------------------------------
    # TP_orbit = 29499.67820209 * u.s # (90 * u.min).to(u.s)
    # m_0 = 5000 * u.kg
    TP_orbit = TP_orbit.to(u.s)
    m_0 = m_0.to(u.kg)
    v_f = np.sqrt(G*M_Earth/to_alt.to(u.m))
    
    I_sp = 4000 * u.second
    T = 400e-3 * u.kg * u.m  / (u.second * u.second) # newtons
    # v_f = 11.2e3 * u.m / u.second
    #---------------------------------------------------

    R_orbit = (TP_orbit*TP_orbit*G*M_Earth / (4*np.pi*np.pi))**(1/3)

    v_0 = np.sqrt(G*M_Earth/R_orbit)

    # v_f = v_0 + delta_v

    v_e = g_0 * I_sp

    p = T / v_e

    m_f = m_0 * np.exp(-1 * (v_f - v_0) / v_e)

    Time_Elapsed = (m_0 - m_f)/p

    a_0 = -1 * p * v_e/ m_0
    a_f = -1 * p * v_e/ m_f

    print(m_0, m_f)
    print(a_0, a_f)
    print(v_0, v_f)
    print(Time_Elapsed, Time_Elapsed.to(u.year))


def new_get_cost(m_0, delta_V):

    #----------------------------------------------------
    g_0 = 9.80665 * u.m / (u.second * u.second)
    G = 6.67408e-11 * u.m**3 * u.kg**(-1) * u.s**(-2)
    M_Earth = 5.972e24 * u.kg
    R_Earth = (6371 * u.km).to(u.m)
    #----------------------------------------------------

    #----------------------------------------------------
    # TP_orbit = 29499.67820209 * u.s # (90 * u.min).to(u.s)
    # m_0 = 5000 * u.kg
    m_0 = m_0.to(u.kg)
    delta_V = delta_V.to(u.m/u.s)
    
    I_sp = 4000 * u.second
    T = 400e-3 * u.kg * u.m  / (u.second * u.second) # newtons
    # v_f = 11.2e3 * u.m / u.second
    #---------------------------------------------------

    v_e = g_0 * I_sp

    p = -1 * T / v_e

    m_f = m_0 * np.exp(-1 * (delta_V) / v_e)

    Time_Elapsed = (m_f - m_0)/p

    a_0 = -1 * p * v_e/ m_0
    a_f = -1 * p * v_e/ m_f

    print(m_0, m_f)
    print(a_0, a_f)
    print(delta_V)
    print(Time_Elapsed, Time_Elapsed.to(u.year))


# new_get_cost(5000*u.kg, 1146.14916511*u.m/u.s)

# new_get_cost(4856.020*u.kg, 888.86125395*u.m/u.s)

# new_get_cost(4747.22*u.kg, 10.80*u.km/u.s)
# new_get_cost(3604.69*u.kg, 7.65*u.km/u.s)

# new_get_cost(2966.00*u.kg, 0.7426167361242552*u.km/u.s)

new_get_cost(2910.38*u.kg, 0.9285591733102699*u.km/u.s)