import pandas as pd
import numpy as np
from astropy import time
from astropy import units as u
from poliastro.twobody import Orbit
from poliastro.bodies import Sun

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

def calc_distance(r1, r2):

    deltax = r2[0].value - r1[0].value
    deltay = r2[1].value - r1[1].value
    deltaz = r2[2].value - r1[2].value

    # deltax *= r1[0].unit
    # deltay *= r1[0].unit
    # deltaz *= r1[0].unit
    
    final_r = np.sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz)

    final_r = final_r * u.km

    return final_r.to(u.au)
