import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.plotting.static import StaticOrbitPlotter
import numpy as np
from astropy import time
import mission_header as mission

Earth_P = mission.Planet("Earth")
UPS_epoch = time.Time("2021-05-01 00:00:00", scale="utc")
Earth_P = Earth_P.propagate(UPS_epoch)

R_orbit = 6652420.04287553 * u.m
probe = Orbit.circular(Earth, alt=R_orbit)
probe = probe.propagate(UPS_epoch)

op = StaticOrbitPlotter()
op.plot(probe)
# plt.show()

