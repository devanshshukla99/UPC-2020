from header import *

"""
Take to GTO
"""

Earth_P = mission.Planet("Earth")
UPS_epoch = time.Time("2021-05-01 00:00:00", scale="utc")
Earth_P = Earth_P.propagate(UPS_epoch)

R_orbit = 6652420.04287553 * u.m
probe = Orbit.circular(Earth, alt=R_orbit)
probe = probe.propagate(UPS_epoch)

mn_0 = Maneuver.hohmann(probe, 32756*u.km)

probe_mn0_a, probe_mn0_f = probe.apply_maneuver(mn_0, intermediate=True)

op = StaticOrbitPlotter()
op.plot(probe)
op.plot(probe_mn0_f)
op.plot(probe_mn0_a)