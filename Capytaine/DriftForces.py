import numpy as np

def drift_forces_delhommeau(data,ckoch,wavamp,numfreq,numdir,numtet):
        
    Fx = np.zeros((numfreq,numdir)).astype(float)
    Fy = np.zeros((numfreq,numdir)).astype(float)
    
    for i in range(numfreq):  
        for j in range(numdir):
            # Pre-calculation of parameters
            # Calculation of deep water wavenumber
            wavnumb0 = (data['omega'][i]**2)/9.81
            # -2* pi * rho * A * w  
            c1 = -2 * np.pi * data['rho']* wavamp * data['omega'][i] 
            # -8* pi * rho * k(k0*h)^2 / (h * [ (k*h)^2 - (k0*h)^2 + k0*h])  
            c2 = -8 * np.pi * data['rho'] * data['wavenumber'][i]*(wavnumb0*data['water_depth'])**2       
            c2 = c2/(data['water_depth']*((data['wavenumber'][i]*data['water_depth'])**2
            - (wavnumb0*data['water_depth'])**2 + wavnumb0*data['water_depth']))

            # ind_theta used for determining the appropriate angle of kochin(w,dir,theta)
            ind_theta = -1
            for k in range(numtet):
                if abs(data["theta"][k]-data["wave_direction"][j]) == np.min(abs(data["theta"]-data["wave_direction"][j])):
                    ind_theta = k
                    break
            if ind_theta != -1:
                print("element index  : ", ind_theta)
            else:
                print("The element not present")

            # The drift forces are a combination of radiation and diffraction kochin functions
            # FORMULA (2.170) in Delhommeau Thesis
            #  c1 * cos(beta) * Im(kochin(w,dir,theta)) + 
            #  c2 * trapezoid(Real(kochin(w,dir,theta))*Imag(Conj(kochin(w,dir,theta)))*cos(theta)::d_theta)  
    
            Fx[i][j] = c1*np.cos(data["wave_direction"][j])*np.imag(ckoch[i][j][ind_theta]) + \
                       c2*np.trapz(np.real(ckoch[i][j][:])*np.imag(np.conj(ckoch[i][j][:]))*np.cos(data["theta"]), x=data["theta"])
            
            Fy[i][j] = c1*np.sin(data["wave_direction"][j])*np.imag(ckoch[i][j][ind_theta]) + \
                       c2*np.trapz(np.real(ckoch[i][j][:])*np.imag(np.conj(ckoch[i][j][:]))*np.sin(data["theta"]), x=data["theta"])
            
    return Fx,Fy
