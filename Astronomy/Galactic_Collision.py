## This program is part of the following publication:
## "W. Schmidt, M. Völschow, Numerical Python in Astronomy and Astrophysics,
## A practical guide of astrophysical problem solving, Springer, 2021"
## See https://link.springer.com/book/10.1007/978-3-030-70347-9


## Galactic collisions

# Module ```galcol``` contains functions for setup, numerical integration, and visualization

# In[ ]:


import galcol


# In[ ]:


dir(galcol)


# #### Example: Whirlpool-like galaxy

# This is an example for a nearly edge-on collision. A smaller intruder galaxy moves under an angle of 45° 
# in the $xy$-plane with an impact velocity of about 130 km/s toward a larger, more massive target galaxy. 
# The impact parameter is 6 kpc in $z$-direction.

# First we define the parameters of the intruder and target galaxies.

# In[ ]:


import galcol
import astropy.units as unit


# In[ ]:


galaxies = {
    'intruder' : galcol.parameters(
        # mass in solar masses
        1e10, 
        # disk radius in kpc
        5, 
        # Cartesian coordinates (x,y,z) of initial position in kpc 
        (25,-25,-5), 
        # x-, y-, z-components of initial velocity in km/s
        (-75,75,0),
        # normal to galactic plane (disk is in xy-plane)
        (0,0,1),
        # number of rings (each ring will be randomly populated with 1000/5 = 200 stars)
        5, 
        # total number of stars
        1000, 
        # softening factor defines inner edge of disk (in units of disk radius)
        0.025),
    'target' : galcol.parameters(5e10, 10, (-5,5,1), (15,-15,0), (1,-1,2**0.5), 10, 4000, 0.025),
}


# In[ ]:


galaxies['intruder']


# In[ ]:


galcol.init_disk(galaxies['intruder'])
galcol.init_disk(galaxies['target'])


# The dictionaries of the two galaxies now contain additional items, 
# particularly the inititial data of the stars.

# In[ ]:


galaxies['intruder']


# Solve equations of motion and visualize data.

# In[ ]:


t, data = galcol.evolve_two_disks(galaxies['target'], galaxies['intruder'], 
                                  N_steps=10000, N_snapshots=500, time_step=0.05*unit.Myr)


# In[ ]:


i = 100
galcol.show_two_disks_3d(data[i,:,:], galaxies['target']['N_stars'], 
                         [-15,15], [-15,15], [-15,15], t[i], name='two_disks')


# In[ ]:


galcol.anim_two_disks_3d(data, galaxies['target']['N_stars'], 
                         [-15,15], [-15,15], [-15,15], t, name='two_disks')


