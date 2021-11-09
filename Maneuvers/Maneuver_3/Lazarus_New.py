
# coding: utf-8

# In[380]:


import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun, Saturn
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.plotting.static import StaticOrbitPlotter
from poliastro.plotting import OrbitPlotter3D
import numpy as np
from astropy import time
from poliastro.twobody import thrust
from poliastro.twobody.propagation import cowell
from poliastro.ephem import Ephem
from poliastro.util import norm, time_range


# In[381]:


date_launch = time.Time("2024-6-26 00:01", scale="utc").tdb # Starting Mission from parking orbit
date_arrival = time.Time("2034-11-15 00:01", scale="utc").tdb


# In[382]:


earth = Ephem.from_body(Earth, time_range(date_launch,end=date_arrival))
saturn = Ephem.from_body(Saturn, time_range(date_launch,end=date_arrival))


# In[383]:


op = StaticOrbitPlotter()
op.plot_body_orbit(Earth, date_launch, label="Earth")


# In[384]:


ss_earth = Orbit.from_ephem(Sun, earth, date_launch)
ss_saturn = Orbit.from_ephem(Sun, saturn, date_arrival)

print(ss_saturn.rv()[0].to(u.km))
print(ss_saturn.rv()[1].to(u.km/u.s))


# In[385]:


op = StaticOrbitPlotter()
op.plot(ss_earth, label="Earth")
op.plot(ss_saturn, label="Saturn")


# In[282]:


man_lambert = Maneuver.lambert(ss_earth, ss_saturn) # Maneuver II


# In[288]:


ss_trans, ss_target = ss_earth.apply_maneuver(man_lambert, intermediate=True)


# In[366]:


op = StaticOrbitPlotter()
op.plot(ss_earth, label="Earth")
op.plot(ss_trans, label="Transfer Orbit")
op.plot(ss_saturn, label="Saturn")


# In[365]:


plotter = OrbitPlotter3D() # StaticOrbitPlotter() #
plotter.set_attractor(Sun)

plotter.plot_body_orbit(Earth, date_launch, label="Earth at launch position")
plotter.plot_body_orbit(Saturn, date_arrival, label="Saturn at arrival position")#, color="C1")
plotter.plot_trajectory(ss_trans.sample(max_anomaly=350*u.deg), label="Tranfer Orbit")#, color="C2")


# In[364]:


#Impulse
man_lambert.impulses[1][0]


# In[362]:


#Impulse
print(man_lambert.impulses[0][1].to(u.km/u.s))
print(norm(man_lambert.impulses[0][1].to(u.km/u.s)))


# In[363]:


print(man_lambert.impulses[1][1].to(u.km/u.s))
print(norm(man_lambert.impulses[1][1].to(u.km/u.s)))


# In[386]:


# Maneuver 3


# In[387]:


ss_trans_in = ss_trans.propagate(date_arrival)


# In[388]:


ss_trans_in_sat = ss_trans_in
ss_trans_in_sa = ss_trans_in_sat.change_attractor(Saturn, force=True)


# In[389]:


ss_ = Orbit.from_ephem(Saturn, saturn, date_arrival)
ss_final_inc = ss_.circular(Saturn, alt=271e5*u.km).propagate(date_arrival)


# In[390]:


op = StaticOrbitPlotter()

op.plot(ss_final_inc, label="Final Orbit")
op.plot(ss_trans_in_sa, label="Tranfer Orbit")
# op.plot(ss_saturn_i, label="Probe around Saturn")


# In[391]:


# ss_trans_in_sa = ss_trans_in_sa.propagate(1*u.day)
# # man_inc_final = Maneuver.lambert(ss_final_inc, ss_trans_in_sa)
man_inc_final = Maneuver.hohmann(ss_trans_in_sa, 271e5*u.km )

ss_trans_final, ss_final = ss_trans_in_sa.apply_maneuver(man_inc_final, intermediate=True)


# In[392]:


op = StaticOrbitPlotter()

op.plot(ss_trans_in_sa, label="Probe Orbit from Maneuver II")
op.plot(ss_trans_final, label="Transfer Orbit")
op.plot(ss_final, label="Final Orbit")


# In[393]:


man_inc_final.impulses


# In[394]:


print(man_inc_final.impulses[1][1].to(u.km/u.s))
print(norm(man_inc_final.impulses[1][1].to(u.km/u.s)))

