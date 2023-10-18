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


# spatial libraries
import numpy as np

class VesselProperties:
    """Mixin class: Something that has vessel properties
    This mixin is updated to better accommodate the ConsumesEnergy mixin

    - type: can contain info on vessel type
    - L: vessel length between perpendiculars (m)
    - B: vessel width (m)
    - Tb: bow draught (m)
    - Ts: stern draught (m)
    - Disp: load displacement (t)
    - C_WP: waterplane coefficient [-] (If none, calculated)
    - C_M: midship section coefficient [-] (If none, calculated)
    - Npro: number of propellers
    - bulbous_bow: False/True (inland ships generally do not have a bulbous_bow, 
    set to False (default). If a ship has a bulbous_bow, set to True)
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
     For a restricted channel: W < Weff; height hT < h0+Dwl; slope Nb >=0
     For a canal: W < Weff; height hT = h0+Dwl; slope Nb>=0)  
    """

    def __init__(
        self,
        type,
        L,
        B,
        Tb,
        Ts,
        Displ,
        C_WP,
        C_M,
        Npro,
        h0,
        hT,
        W,
        Nb,
        Dwl = 0.0,
        rho = 1.0,
        bulbous_bow = False,
        transom_stern = False,
        safety_margin=None,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)

        """Initialization
        """
        self.type = type
        self.B = B
        self.L = L
        self.Tb = Tb
        self.Ts = Ts
        self.Displ = Displ
        self.C_WP = C_WP
        self.C_M = C_M
        self.Npro = Npro
        self.bulbous_bow = bulbous_bow
        self.transom_stern = transom_stern

        self.Dwl = Dwl       
        self.rho = rho
        self.h0 = h0
        self.hT = hT
        self.W = W
        self.Nb = Nb
        
        # alternative  options
        self.safety_margin = safety_margin
        

    def calculation_init (self):
        """Calculation of initial values:
        - Tm: mean draught (m)
        - Ukc: underkeel clearance = (h0+Dwl)/Tm
        - Weff: channel effective width (m)
        - C_B: block coefficient
        - As: midship section (m2)
        - Ach: channel section (m2)
        """
        self.Tm = (self.Tb+self.Ts)/2
        self.Ukc = (self.h0+self.Dwl)/self.Tm
        self.C_B = self.Displ/(self.L*self.B*self.Tm*self.rho)
        self.Disp = self.Displ/self.rho
        self.Weff = 7.04*self.B/(self.C_B**0.85)
        if self.C_WP is None:
            self.C_WP = (1 + 2 * self.C_B) / 3  # Waterplane coefficient            
        if self.C_M is None:
            self.C_M = 1.006 - 0.0056 * self.C_B ** (-3.56)  # Midship section coefficient
        self.As = self.C_M*self.B*self.Tm
        """ hT should be less or equal to h0+Dwl"""
        if (self.h0+self.Dwl) < self.hT:
            self.hT = self.h0+self.Dwl
        if self.W <= self.Weff:
            assert self.hT > 0 and self.hT <= (self.h0+self.Dwl), f"W less than Weff,then hT should be positive and less or equal h0+Dwl: {self.hT}"
            assert self.Nb >= 0, f"Nb should be positive or 0: {self.Nb}"
            self.Ach = (self.W+self.Nb*(self.h0+self.Dwl))*(self.h0+self.Dwl)
        else:
            assert self.hT == 0, f"hT should be 0: {self.hT}"
            assert self.Nb == 0, f"Nb should be 0: {self.Nb}"
            self.Ach = self.Weff*(self.h0+self.Dwl)
        assert self.Npro > 0, f"Npro should be greater than 0: {self.Npro}" 
        """Calculation of critical speed Vcr:
        - Kch: coefficient (unrestricted channel)
        - Kc: coefficient (restricted channel)
        - hm: mean water depth (rectangular section)
        - hmT: mean water depth (restricted channel)
        """
        self.Kch = 0.58*((self.h0+self.Dwl)*self.L/self.B/self.Tm)**0.125
        self.Kc = (2*np.cos((np.pi+np.arccos(1-self.As/self.Ach))/3))**1.5
        if self.W <= self.Weff:
            self.hm = self.Ach/(self.W+2*self.Nb*(self.h0+self.Dwl))
        else:
            self.hm = self.Ach/(self.Weff+2*self.Nb*(self.h0+self.Dwl))
        self.hmT = (self.h0+self.Dwl)-self.hT*(1-self.hm/(self.h0+self.Dwl))
        if self.hT == 0:
            self.Vcr = self.Kch*np.sqrt(9.81*(self.h0+self.Dwl))
        elif self.hT < (self.h0+self.Dwl):
            self.Vcr = (self.Kch*(1-self.hT/(self.h0+self.Dwl))+self.Kc*self.hT/(self.h0+self.Dwl))*np.sqrt(9.81*self.hmT)
        elif self.hT == (self.h0+self.Dwl):
            self.Vcr = self.Kc*np.sqrt(9.81*self.hmT)
        """Calculation of limit speed Vlim:
        - Vlim: limit speed = 0.9*Vcr
        """
        #self.Vlim = 0.9*self.Vcr
        self.Vlim = self.Vcr
        
        return self.Weff, self.Vlim 

    def calculation_Hooft(self,v):
        """Squat calculation for unrestricted channel (Hooft, 1974):
        - Fnh: Froude number
        - Shb: squat at the bow (m)
        """
        self.Fnh = v/np.sqrt(9.81*(self.h0+self.Dwl))
        if self.hT == 0:
            self.Shb = 2*self.C_B*self.B*self.Tm*(self.Fnh**2)/self.L/np.sqrt(1-self.Fnh**2)
            if self.Ukc < 1.2 and self.C_B >= 0.8: self.Shb = 0
        else:
            self.Shb = 0
        
        return self.Shb
    
    def calculation_Romisch(self,v):
        """Squat calculation for unrestricted channel and canal (Römisch, 1989):
        - Kdt: coefficient
        - Cf: coeficient (different from bow and stern)
        - Cv: coefficient
        - Srb: squat at the bow (m)
        - Srs: squat at the stern (m)
        - Vnoeu: speed in knots
        """
        self.Kdt = 0.155*np.sqrt((self.h0+self.Dwl)/self.Tm)
        self.Cf = (10*self.B*self.C_B/self.L)**2
        self.Cv = 8*np.power(v/self.Vcr,2)*(0.0625+np.power(v/self.Vcr-0.5,4))
        if self.hT == 0:
            self.Srb = self.Cv*self.Cf*self.Kdt*self.Tm
            if self.Ukc < 1.2 and self.C_B < 0.8: self.Srb = 0
            self.Srs = self.Srb/self.Cf
        elif self.hT < (self.h0+self.Dwl):
            self.Srb = 0
            self.Srs = 0
        elif self.hT == (self.h0+self.Dwl):
            self.Srb = self.Cv*self.Cf*self.Kdt*self.Tm
            Vnoeu = v*3600/1852
            if self.Ukc < 1.2 and self.C_B > 0.8 and Vnoeu < 7: self.Srb = 0
            self.Srs = self.Srb/self.Cf

        #print(self.B,self.C_B,self.Cf,self.Srb)
        
        return self.Srb, self.Srs

    def calculation_Ankudinov(self,v):
        """Squat calculation for restricted channel and canal (Ankudinov, 2000):
        - Kps, Kpt: propeller coefficients
        - Phu, Pfnh, Pht, Sh, Pch1: coefficients
        - Kbt: forward bulb coefficient
        - Ktrt: transom stern coefficient
        - Kt1t, Ktr, Phtm, Pch2: coefficients 
        - Sab: squat at the bow (m)
        - Sas: squat at the stern (m)
        - Trim: ship trim (m)
        - Vnoeu: speed in knots
        """
        if self.Npro == 1:
            self.Kps = 0.15
        else:
            self.Kps = 0.13
        self.Phu = 1.7*self.C_B*(self.B*self.Tm/self.L**2)+0.004*self.C_B**2
        self.Fnh = v/np.sqrt(9.81*(self.h0+self.Dwl))
        self.Pfnh = np.power(self.Fnh,1.8+0.4*self.Fnh)
        self.Pht = 1+0.35/self.Ukc**2
        self.Sh = self.C_B*self.Tm*self.hT*self.As/self.Ach/(self.h0+self.Dwl)**2
        if self.hT == 0:
            self.Pch1 = 1
        else:
            self.Pch1 = 1+10*self.Sh-1.5*(1+self.Sh)*np.sqrt(self.Sh)
        self.Sab = self.L*(1+self.Kps)*self.Phu*self.Pfnh*self.Pht*self.Pch1

        #print(self.L,self.Kps,self.Phu,self.Pfnh,self.Pht,self.Pch1,self.Sab)        

        if self.Npro == 1:
            self.Kpt = 0.15
        else:
            self.Kpt = 0.20
        if self.bulbous_bow:
            self.Kbt = 0.1
        else:
            self.Kbt = 0
        if self.transom_stern:
            self.Ktrt = 0.04
        else:
            self.Ktrt = 0        
        self.Kt1t = (self.Ts - self.Tb)/(self.Ts+self.Tb)
        self.Ktr = np.power(self.C_B,2+0.8*self.Pch1/self.C_B)-(0.15*self.Kps+self.Kpt)-(self.Kbt+self.Ktrt+self.Kt1t)
        self.Phtm = 1-np.exp(2.5*(1-self.Ukc)/self.Fnh)
        if self.hT == 0:
            self.Pch2 = 1
        else:
            self.Pch2 = 1-5*self.Sh
        self.Trim = -1.7*self.L*self.Phu*self.Pfnh*self.Phtm*self.Ktr*self.Pch2        

        #print(self.L,self.Phu,self.Pfnh,self.Phtm,self.Ktr,self.Pch2,self.Trim)        
        
        if self.hT == 0:
            self.Sab = 0
            self.Sas = 0
        elif self.hT < (self.h0+self.Dwl):
            self.Sas = self.Sab+0.5*self.Trim
            self.Sab = self.Sab-0.5*self.Trim
        elif self.hT == (self.h0+self.Dwl):
            self.Sas = self.Sab+0.5*self.Trim
            self.Sab = self.Sab-0.5*self.Trim
            Vnoeu = v*3600/1852
            if self.Ukc < 1.2 and self.C_B > 0.8 and Vnoeu > 7:
                self.Sab = 0
                self.Sas = 0
        
        
        return self.Sab, self.Sas

    def calculate_squat(self, v):
        """Squat for unrestricted channel, restricted channel or canal:
        - Squat_v: squat value for a ship speed 
        """
        
        self.calculation_init()
        self.Shb = self.calculation_Hooft(v)
        self.Srb,self.Srs = self.calculation_Romisch(v)
        self.Sab,self.Sas = self.calculation_Ankudinov(v)
        
        # The squat is determined according to the type of channel or canal
        if self.hT == 0:
            list1 = [self.Shb,self.Srb,self.Srs]
            self.Squat_v = np.max(list1)
        elif self.hT < (self.h0+self.Dwl):
            list1 = [self.Sab,self.Sas]
            self.Squat_v = np.max(list1)
        elif self.hT == (self.h0+self.Dwl):
            list1 = [self.Sab,self.Sas,self.Srb,self.Srs]
            self.Squat_v = np.max(list1)
            #print("The type is : ",type(self.Squat_v))        
        return self.Squat_v

        
class ExtraMetadata:
    """store all leftover keyword arguments as metadata property (use as last mixin)"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args)
        # store all other properties as metadata
        self.metadata = kwargs
