import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.plotting.static import StaticOrbitPlotter
import numpy as np
from astropy import time
import mission_header as mission
from poliastro.twobody import thrust
from poliastro.twobody.propagation import cowell
import tqdm

Earth_P = mission.Planet("Earth")
UPS_epoch = time.Time("2021-05-01 00:00:00", scale="utc")
Earth_P = Earth_P.propagate(UPS_epoch)
                                                                            
R_orbit = 6652420.04287553 * u.m
probe = Orbit.circular(Earth, alt=R_orbit)
probe = probe.propagate(UPS_epoch)



def updatePosition(rivi, delta):
    # Return tuple with updated position
    # return tuple(i + delta * vi for i, vi in zip(self.position, self.velocity)) # for i in x, y

    ri = rivi[0]
    vi = rivi[0]

    rf = [0.0, 0.0, 0.0]

    rf[0] = ri[0].value + delta * vi[0].value
    rf[1] = ri[1].value + delta * vi[1].value
    rf[2] = ri[2].value + delta * vi[2].value

    return rf * ri.unit

def updateVelocity(rivi, delta):
    # Return tuiple with updated velocity, using leap-frog integration scheme
    # return tuple(i - (delta * G * M * ri) / self.distanceToSun()**3 for i, ri in zip(self.velocity, self.position))

    ri = rivi[0]
    vi = rivi[0]
    
    r_mag = np.sqrt(ri[0]**2 + ri[1]**2 + ri[2]**2)

    vi_f = [0.0, 0.0, 0.0]

    vi_f[0] = vi[0] - delta * G * M * ri[0]/r_mag**3
    vi_f[1] = vi[1] - delta * G * M * ri[1]/r_mag**3
    vi_f[2] = vi[2] - delta * G * M * ri[2]/r_mag**3

    return vi_f * vi.unit

def updaterv(rivi, delta):

    M_Earth = 5.972e24 * u.kg

    ri = rivi[0]
    vi = rivi[0]

    rf = [0.0, 0.0, 0.0]

    delta = delta.to(u.s)

    rf[0] = ri[0].value + delta.value * vi[0].value
    rf[1] = ri[1].value + delta.value * vi[1].value
    rf[2] = ri[2].value + delta.value * vi[2].value

    r_mag = np.sqrt(ri[0]**2 + ri[1]**2 + ri[2]**2).to(u.m)

    vi_f = [0.0, 0.0, 0.0]

    vi_f[0] = vi[0] - (delta * G * M_Earth * ri[0]/r_mag**3)
    vi_f[1] = vi[1] - (delta * G * M_Earth * ri[1]/r_mag**3)
    vi_f[2] = vi[2] - (delta * G * M_Earth * ri[2]/r_mag**3)

    return rf * ri.unit, vi_f



class Probe():

    def __init__(self):
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

        self.g_0 = 9.80665 * u.m / (u.second * u.second)
        self.G = 6.67408e-11 * u.m**3 * u.kg**(-1) * u.s**(-2)
        self.M_Earth = 5.972e24 * u.kg
        self.R_Earth = (6371 * u.km).to(u.m)

        self.TP_orbit_Initial = (90 * u.min).to(u.s)
        self.TP_orbit = self.TP_orbit_Initial

        self.m_0 = 5000 * u.kg

        self.I_sp = 4000 * u.second
        self.Thrust = 400e-3 * u.kg * u.m  / (u.second * u.second) # newtons

        self.R_orbit = ( self.TP_orbit_Initial*self.TP_orbit_Initial*self.G*self.M_Earth / (4*np.pi*np.pi) )**(1/3)
        self.v_0 = np.sqrt(self.G*self.M_Earth/self.R_orbit)

        self.probe = Orbit.circular(Earth, alt=R_orbit)
        self.UPS_epoch = time.Time("2021-05-01 00:00:00", scale="utc")
        self.probe = probe.propagate(self.UPS_epoch)

        self.v_e = self.g_0 * self.I_sp

        self.rate_of_mass = -1 * self.Thrust / self.v_e

        self.m_inst = self.m_0         # Inst t

    # def give_thrust(self, delta_t):f

    #     delta_t = delta_t.to(u.s)

    #     self.v_f = [0.0, 0.0, 0.0]

    #     self.r_i = self.probe.rv()[0]
    #     self.v_i = self.probe.rv()[1]

    #     self.m_inst = self.m_inst + self.rate_of_mass * delta_t

    #     self.v_f[0] = self.v_i[0] + self.v_e * np.log(self.m_0/self.m_inst)
    #     self.v_f[1] = self.v_i[1] + self.v_e * np.log(self.m_0/self.m_inst)
    #     self.v_f[2] = self.v_i[2] + self.v_e * np.log(self.m_0/self.m_inst)
    
    def give_thrust(self, time_of_burn):
        
        #----------------------------------------------------#
        self.m_inst = (self.m_inst).to(u.kg)
        # self.TP_orbit = self.TP_orbit.to(u.s)
        time_of_burn = (time_of_burn).to(u.s)

        # self.R_orbit = ( (self.TP_orbit*self.TP_orbit*self.G*self.M_Earth) / (4*np.pi*np.pi) )**(1/3)
        # self.v_orbit = np.sqrt(self.G*self.M_Earth/self.R_orbit)

        self.v_orbit = self.norm_(self.probe.v)

        self.m_inst = self.m_inst + time_of_burn*self.rate_of_mass
        self.v_f = self.v_orbit + self.v_e * np.log(self.m_0/self.m_inst)

    def norm_(self, quantity):
        return np.sqrt(quantity[0]*quantity[0] + quantity[1]*quantity[1] + quantity[2]*quantity[2])

    def update_rv_after_1_rev(self, time_of_burn):
        
        time_of_burn = (time_of_burn).to(u.s)

        self.r_i = self.probe.rv()[0]
        self.v_i = self.probe.rv()[1]

        delta_change_in_v = (self.v_f**2 - self.norm_(self.v_i)**2) / ( 2* (self.v_i[0] + self.v_i[1]) ) 

        self.v_i[0]  += delta_change_in_v
        self.v_i[1]  += delta_change_in_v

        self.probe.v = self.v_i

        self.r_f = [0.0, 0.0, 0.0]

        __lx__ = time_of_burn.value/1e5

        for __l__ in tqdm.tqdm(range(0, 100000)):
            self.r_i[0] = (self.r_i[0].to(u.m).value + __lx__ * self.v_i[0].to(u.m/u.s).value) * u.m
            self.r_i[1] = (self.r_i[1].to(u.m).value + __lx__ * self.v_i[1].to(u.m/u.s).value) * u.m
            self.r_i[2] = (self.r_i[2].to(u.m).value + __lx__ * self.v_i[2].to(u.m/u.s).value) * u.m

        self.r_f[0] = self.r_i[0].value
        self.r_f[1] = self.r_i[1].value
        self.r_f[2] = self.r_i[2].value

        self.probe.r = (self.r_f * u.m).to(u.km)

        # print(self.r_f)

    def old_update_rv_after_1_rev(self, time_of_burn):

        time_of_burn = (time_of_burn).to(u.s)

        self.r_i = self.probe.rv()[0]
        self.v_i = self.probe.rv()[1]

        delta_change_in_v = (self.v_f**2 - self.norm_(self.v_i)**2) / ( 2* (self.v_i[0] + self.v_i[1]) ) 

        self.v_i[0]  += delta_change_in_v
        self.v_i[1]  += delta_change_in_v

        self.probe.v = self.v_i

        self.r_f = [0.0, 0.0, 0.0]

        self.r_f[0] = (self.r_i[0].to(u.m).value + time_of_burn.value * self.v_i[0].to(u.m/u.s).value)
        self.r_f[1] = (self.r_i[1].to(u.m).value + time_of_burn.value * self.v_i[1].to(u.m/u.s).value)
        self.r_f[2] = (self.r_i[2].to(u.m).value + time_of_burn.value * self.v_i[2].to(u.m/u.s).value)

        self.probe.r = (self.r_f*u.m).to(u.km)

        print(self.r_f)

    def infinitisimal_thrust(self, time_of_burn):

        self.m_inst = (self.m_inst).to(u.kg)

        time_of_burn = (time_of_burn).to(u.s)

        self.v_orbit = self.norm_(self.probe.v)

        self.r = self.probe.rv()[0]
        self.v = self.probe.rv()[1]

        self.r_f = [0.0, 0.0, 0.0]

        for ti in tqdm.tqdm(range(0, int(time_of_burn.value))):

            self.m_inst = (self.m_inst).to(u.kg) + ti * (self.rate_of_mass).to(u.kg/u.s) * u.s
            self.v_f = (self.norm_(self.v).to(u.m/u.s)) + (self.v_e.to(u.m/u.s)) * np.log( (self.m_0).to(u.kg)/(self.m_inst).to(u.kg) )

            self.acceleration = -1 * (self.rate_of_mass.to(u.kg/u.s)) * (self.v_e.to(u.m/u.s)) / (self.m_inst.to(u.kg))
            
            self.delta_change_in_v = ( self.v_f.to(u.m/u.s)**2 - (self.norm_(self.probe.v).to(u.m/u.s))**2 ) / ( 2* ( self.probe.v[0].to(u.m/u.s) + self.probe.v[1].to(u.m/u.s) ) )

            self.v[0]  += self.delta_change_in_v.to(u.km/u.s)
            self.v[1]  += self.delta_change_in_v.to(u.km/u.s)

            self.probe.v = self.v

            self.r_f[0] = (self.r[0].to(u.m).value + ti*self.v[0].to(u.m/u.s).value + 0.5 * ti * ti * self.acceleration.to(u.m/u.s/u.s).value)
            self.r_f[1] = (self.r[1].to(u.m).value + ti*self.v[1].to(u.m/u.s).value + 0.5 * ti * ti * self.acceleration.to(u.m/u.s/u.s).value)
            self.r_f[2] = (self.r[2].to(u.m).value + ti*self.v[2].to(u.m/u.s).value)

            self.r = [0.0, 0.0, 0.0]

            self.r[0] = self.r_f[0]
            self.r[1] = self.r_f[1]
            self.r[2] = self.r_f[2]

            self.r = self.r * u.m

            self.probe.r = self.r.to(u.km)


        
    pass