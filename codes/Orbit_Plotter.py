import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.plotting.static import StaticOrbitPlotter
import pandas as pd

"""
  Symbol meaning [1 au= 149597870.700 km, 1 day= 86400.0 s]:

    JDTDB    Julian Day Number, Barycentric Dynamical Time
      EC     Eccentricity, e
      QR     Periapsis distance, q (au)
      IN     Inclination w.r.t X-Y plane, i (degrees)
      OM     Longitude of Ascending Node, OMEGA, (degrees)
      W      Argument of Perifocus, w (degrees)
      Tp     Time of periapsis (Julian Day Number)
      N      Mean motion, n (degrees/day)
      MA     Mean anomaly, M (degrees)
      TA     True anomaly, nu (degrees)
      A      Semi-major axis, a (au)
      AD     Apoapsis distance (au)
      PR     Sidereal orbit period (day)
"""

def read_jpl_file(file_):
        
    df = pd.read_csv(file_ + ".jpl", delimiter="=", index_col=0, header=None) 
    data = df.to_dict()[1]

    data['EC'] = float(data['EC' ].strip(" "))
    data['QR'] = float(data['QR' ].strip(" "))
    data['IN'] = float(data['IN' ].strip(" "))
    data['OM'] = float(data['OM' ].strip(" "))
    data['W '] = float(data['W ' ].strip(" "))
    data['Tp'] = float(data['Tp' ].strip(" "))
    data['N '] = float(data['N ' ].strip(" "))
    data['MA'] = float(data['MA' ].strip(" "))
    data['TA'] = float(data['TA' ].strip(" "))
    data['A '] = float(data['A ' ].strip(" "))
    data['AD'] = float(data['AD' ].strip(" "))
    data['PR'] = float(data['PR'].strip(" "))

    return data

def Planet(planet_name):

    data = read_jpl_file(planet_name)

    return Orbit.from_classical(
        attractor=Sun,
        a = data["A "] * u.AU,
        ecc = data["EC"] * u.one,
        inc = data["IN"] * u.deg,
        raan = data["OM"] * u.deg,
        argp = data["W "] * u.deg,
        nu = data["TA"] * u.deg
    )

Earth = Planet("Earth")
Mars = Planet("Mars")
Jupiter = Planet("Jupiter")
Saturn = Planet("Saturn")

op = StaticOrbitPlotter()
op.plot(Earth)
op.plot(Mars)
op.plot(Jupiter)
op.plot(Saturn)

# orb = Orbit.
# op=StaticOrbitPlotter()

# # Data for Mars at J2000 from JPL HORIZONS
# a = 1.523679 * u.AU         # SemiMajor Axis
# ecc = 0.093315 * u.one      # Eccentricity
# inc = 1.85 * u.deg          # Inclination
# raan = 49.562 * u.deg       # RA of Ascendng node
# argp = 286.537 * u.deg      # Argument of pericenter
# nu = 23.33 * u.deg          # True Anamoly

# orb = Orbit.from_classical(Sun, a, ecc, inc, raan, argp, nu)

# from poliastro.ephem import Ephem

# from astropy import time
# epoch = time.Time("2020-04-29 10:43")  # UTC by default

# for i in range(10, 20):

#     epoch = time.Time("20" + str(i) + "-04-29 10:43")  # UTC by default
#     op.plot(Mars.propogate(epoch))
