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

Earth_P = mission.Planet("Earth")
UPS_epoch = time.Time("2021-05-01 00:00:00", scale="utc")
Earth_P = Earth_P.propagate(UPS_epoch)
                                                                            
R_orbit = 6652420.04287553 * u.m
probe = Orbit.circular(Earth, alt=R_orbit)
probe = probe.propagate(UPS_epoch)

a_d, _, t_f = thrust.change_ecc_quasioptimal(probe, 0.50, 1e-5 * u.m/u.s/u.s)  

probe_f = probe.propagate(t_f*u.s, method=cowell, ad=a_d)


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



