import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun, Saturn
from poliastro.twobody import Orbit
# from poliastro.maneuver import Maneuver
from poliastro.plotting.static import StaticOrbitPlotter
from poliastro.ephem import Ephem
# import pandas as pd
import numpy as np
from astropy import time

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

    deltax = r2[0].to(u.m) - r1[0].to(u.m)
    deltay = r2[1].to(u.m) - r1[1].to(u.m)
    deltaz = r2[2].to(u.m) - r1[2].to(u.m)

    final_r = np.sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz)

    return final_r.to(u.au)

# Earth = Planet("Earth")
# Mars = Planet("Mars")
# Jupiter = Planet("Jupiter")
# Saturn = Planet("Saturn")

UPS_epoch = time.Time("2034-05-01 00:00:00", scale="utc").tdb

# Earth = Earth.propagate(UPS_epoch)
# Mars = Mars.propagate(UPS_epoch)
# Jupiter = Jupiter.propagate(UPS_epoch)
# Saturn = Saturn.propagate(UPS_epoch)

# earth = Orbit.from_body_ephem(Earth, epoch=UPS_epoch) 
# saturn = Orbit.from_body_ephem(Saturn, epoch=UPS_epoch) 

e_ = Ephem.from_body(attractor=Sun, body=Earth, epochs=UPS_epoch)
s_ = Ephem.from_body(attractor=Sun, body=Saturn, epochs=UPS_epoch)

earth = Orbit.from_ephem(ephem=e_, attractor=Sun, epoch=UPS_epoch)
saturn = Orbit.from_ephem(ephem=s_, attractor=Sun, epoch=UPS_epoch)

integration = 2 * u.day
total_calc = 1*365 * u.day 
dse = []
epochs = []

iterations = total_calc/integration

print("-----------------------------------------------\nInitial\n-----------------------------------------------")
print("Earth", earth.rv())
print("Saturn", saturn.rv())
d_ = calc_distance(earth.rv()[0], saturn.rv()[0])
minimum_d = d_.value

print("DSE", )
print("---------------------------------------------------------------------------------------------------------")

propagate_Earth = earth
propagate_Saturn = saturn

for i in range(0, int(iterations)):

    propagate_Earth = propagate_Earth.propagate(10*u.day)
    propagate_Saturn = propagate_Saturn.propagate(10*u.day)

    d_ = calc_distance(propagate_Earth.rv()[0], propagate_Saturn.rv()[0])
    dse.append(d_.value)

    epochs.append(propagate_Earth.epoch.strftime("%Y-%m-%d %H:%M"))

    if(minimum_d>d_.value):
        print("$ DSE -----> " + str(d_), "\n @ Epoch ---> " + propagate_Earth.epoch.strftime("%Y-%m-%d %H:%M"))


dse = np.array(dse)
np.save("DSE_", dse)
np.save("epochs", epochs)

op = StaticOrbitPlotter()
op.plot(Earth)
op.plot(Saturn)
op.plot(propagate_Earth)
op.plot(propagate_Saturn)
plt.show()
