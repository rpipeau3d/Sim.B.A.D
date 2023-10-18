# -*- coding: utf-8 -*-
"""
This program is free software:
you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the program.
If not, see <https://www.gnu.org/licenses/>.

"""

import Modules.vessel_channel as vech
import numpy as np
import matplotlib.pyplot as plt
import Modules.energy as energy


"""Vessel and channel properties

    - type: can contain info on vessel type
    - L: vessel length between perpendiculars (m)
    - B: vessel width (m)
    - Tb: bow draught (m)
    - Ts: stern draught (m)
    - Displ: load displacement (t)
    - C_WP: waterplane coefficient (If none, calculated)
    - C_M: midship section coefficient (If none, calculated)
    - Npro: number of propellers
    - bulbous_bow: False/True
    - transom_stern: False/True
    
    - rho: density of the surrounding water (t/m^3)
    - Dwl: design water level (m)
    - h0: water depth (m)
    - hT: height of trench (m)
    - W: channel width (m)
    - Nb: inverse bank slope (m/m)
    - safety_margin : the water area above the waterway bed reserved to prevent 
      ship grounding due to ship squatting during sailing,the value of safety 
      margin depends on waterway bed material and ship types. 
      For tanker vessel with rocky bed the safety margin is recommended 
      as 0.3 m based on Van Dorsser et al. The value setting for safety margin
      depends on the risk attitude of the ship captain and shipping companies.

     Note:
     For an unrestricted channel: height hT = 0; slope Nb = 0; width W >= Weff
     For a restricted channel: W < Weff; height 0< hT < h0+Dwl; slope Nb >=0
     For a canal: W < Weff; height hT = h0+Dwl; slope Nb>=0)  

"""
# Make your preferred class out of available mix-ins.
TransportResource = type(
    "Vessel",
    (
        vech.VesselProperties,  # needed to add vessel properties
        energy.ConsumesEnergy,
    ),
    {},
)  # needed to calculate resistances

data_vessel = {"type":None,
               "L" : 205,
               "B" : 32,
               "Tb" : 10, 
               "Ts" : 10, 
               "Displ" : 37500,      
               "C_WP" : 0.75,
               "C_M" : 0.98,
               "Npro" : 1, # number of propellers
               "bulbous_bow" : True, # if a vessel has no bulbous_bow, set to False; otherwise set to True.
               "transom_stern" : True, # if a vessel has no transom stern, set to False; otherwise set to True.
               "rho" : 1.0,
               "Dwl" : 0.0,
               "h0" : 30.0,
               "hT" : 0.0,
               "W" : 400, 
               "Nb" : 0,
               "safety_margin" : 0.2, # for tanker vessel with sandy bed the safety margin is recommended as 0.2 m
              }             

"""Something that consumes energy.

    - P_installed: installed engine power [kW]
    - L_w: weight class of the ship (depending on carrying capacity) 
    (classes: L1 (=1), L2 (=2), L3 (=3))
    - C_year: construction year of the engine [y]
    - P_hotel: power for systems on board [kW]
      (If none, calculated with P_hotel_perc)
    - P_hotel_perc: percentage of P_installed for P_hotel; default value:0.05)
    - nu: kinematic viscosity [m^2/s] (default value:1 * 10 ** (-6) )
    - g: gravitational accelleration [m/s^2] (default value: 9.81)
    - eta_o: open water efficiency of propeller [-] (default value: 0.4)
    - eta_r: relative rotative efficiency [-] (default value:1.O)
    - eta_t: transmission efficiency [-] (default value: 0.98)
    - eta_g: gearing efficiency [-] (default value: 0.96)
    - c_stern: determines shape of the afterbody [-] (default value: 0.0)
    - C_BB: breadth coefficient of bulbous_bow, set to 0.2 according to 
    the paper of Kracht (1970), https://doi.org/10.5957/jsr.1970.14.1.1
    - one_k2: appendage resistance factor (1+k2) [-] (default value: 2.5)
    - Typvessel: type of vessel for calculation of propulsion power (default value: 'Inland')
      'Inland' for riverboat,
      'Tanker' for tankers and bulkcarriers
      'Container' for containerships
      'RoRo' for RoRo ships
    - S_APP1: wetted area of appendages in percentage of hull wetted area (default value: 0.05)
    - h_B1: Position of the centre of the bulb transverse area in percentage of mean draught (default value: 0.2)
    - A_T1: Transverse area of the transom in percentage of B*T (default value: 0.2)
    - D_s: Propeller diameter (If none, calculated)
"""

data_moteur = {"P_installed": 32700.0,
               "L_w": 3.0 ,
               "C_year":2020,
               "c_stern":10,
               "C_BB":0.0638,
               "one_k2":1.5,
               "S_APP1":0.0065,
               "h_B1":0.4 ,
               "A_T1":0.05 ,
               "Typvessel":'Tanker',
               "D_s": 8,
              }             


bateau = TransportResource(**data_vessel,**data_moteur)
Weff, Vlim = bateau.calculation_init()
print("Effective width Weff: {:.2f} m".format(Weff))
print("Limit speed: {:.2f} m/s".format(Vlim))
Squat_v = bateau.calculate_squat(Vlim)
print("Squat: {:.2f} m".format(Squat_v))

# calculate available underkeel clearance (vessel in rest)
z_given = (bateau.h0+bateau.Dwl) - bateau.Tm
Vmax = 20.0
velocity = np.linspace(0.01,Vmax, 1000)
Squat_v = np.zeros(1000)
#print("The type is : ",type(velocity))
# calculate sinkage
Squat_v[0] = bateau.calculate_squat(velocity[0])
# Initialisation
i = 0
# compute difference between the sinkage and the space available for sinkage
diff = z_given - Squat_v[0] - bateau.safety_margin
v = velocity[0]
# The value of velocity within the bound (0,20) is the velocity we find 
# where the diff reach a minimum (zero).
while diff >=0 and v <= Vlim:        
    i+=1
    v = velocity[i]
    Squat_v[i] = bateau.calculate_squat(v)        
    diff = z_given - Squat_v[i] - bateau.safety_margin
    #print(Vessel.hT,velocity[i],Squat_v)
imax = i
grounding_v = velocity[imax]
if grounding_v > Vlim:
    print("Grounding velocity greater than Limit speed")
else:
    print("Grounding velocity: {:.2f} m/s".format(grounding_v))
    print("Squat: {:.2f} m".format(Squat_v[imax]))

# Recalculate squat with min between limit speed and grounding velocity
# Calculate:
# - Total resistance R_tot (kN)
# - Power required to propellers: P_propulsion (kW)
# - Total power required: P_propulsion + P_hotel (kW)
# - Actual total power installed: P_installed (kW)

R_tot = np.zeros(1000)
P_propulsion = np.zeros(1000)
P_tot = np.zeros(1000)
P_installed = np.zeros(1000)
list1 = [Vlim,grounding_v]
Vmax = np.min(list1)
velocity = np.linspace(0.01,Vmax, 1000)
for i in range(0,1000):
    v = velocity[i]
    Squat_v[i] = bateau.calculate_squat(v)
    h_0 = (bateau.h0+bateau.Dwl) - bateau.calculate_squat(v)
    D = h_0 - bateau.Tm
    R_tot[i] = bateau.calculate_total_resistance(v,h_0)
    P_propulsion[i], P_tot[i], P_installed[i] = bateau.calculate_total_power_required(v=v, h_0=h_0)
    bateau.calculate_emission_factors_total(v, h_0)
    bateau.calculate_SFC_final(v, h_0)
  
#Plot of squat
plt.interactive(True)
plt.figure()
plt.plot(velocity[:],Squat_v[:],'r-',label='Squat')
plt.xlabel('Velocity(m/s)')
plt.ylabel('Squat(m)')
plt.grid(True)
plt.legend()
plt.title('Cargo: L$_s$ = 205 m, B$_s$ = 32 m, T$_s$ = 10 m')
plt.savefig('Fig_1.pdf',bbox_inches = 'tight', dpi=600, format='pdf')

#Plot of Total resistance (kN)
plt.interactive(True)
plt.figure()
plt.plot(velocity[:],R_tot[:],'r-',label='Total resistance')
plt.xlabel('Velocity(m/s)')
plt.ylabel('Total resistance(kN)')
plt.grid(True)
plt.legend()
plt.title('Cargo: L$_s$ = 205 m, B$_s$ = 32 m, T$_s$ = 10 m')
plt.savefig('Fig_2.pdf',bbox_inches = 'tight', dpi=600, format='pdf')

#Plot of P_propulsion, P_tot, P_given (kW)
plt.interactive(True)
plt.figure()
plt.plot(velocity[:],P_propulsion[:],'r--',label='P_propulsion')
plt.plot(velocity[:],P_tot[:],'b--',label='P_tot')
plt.plot(velocity[:],P_installed[:],'g-',label='P_installed')
plt.plot(25*1.852/3.6,32621*0.99*0.9931,'ro',label='P installed Holtrop_Mennen')
plt.xlabel('Velocity(m/s)')
plt.ylabel('Power (kW)')
plt.grid(True)
plt.legend()
plt.title('Cargo: L$_s$ = 205 m, B$_s$ = 32 m, T$_s$ = 10 m')
plt.savefig('Fig_3.pdf',bbox_inches = 'tight', dpi=600, format='pdf')
