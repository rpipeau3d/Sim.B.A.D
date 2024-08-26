import capytaine as cpt; print(cpt.__version__)
import numpy as np
import xarray as xr
from DriftForces import drift_forces_delhommeau as DF

cpt.set_logging('INFO')

# Generate the mesh of ship; Keep only the immersed part of the mesh
half_mesh = cpt.load_mesh('Cont_3.stl')
mesh = cpt.ReflectionSymmetricMesh(half_mesh, cpt.xOz_Plane)
mesh = mesh.immersed_part(water_depth=13.10)
# lid_mesh = mesh. generate_lid(z=-0.1)
cog = np.array((180.8, 0, 11.2))
body = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center = cog),
        center_of_mass= cog,
        # lid_mesh = lid_mesh
        )
body.show()
print(body.dofs.keys())
# dict_keys(['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'])

# Hydrostatic calculations
hydrostatics = body.compute_hydrostatics(rho=1025.0)
print("Volume:", hydrostatics["disp_volume"])
print("Center of buoyancy:", hydrostatics["center_of_buoyancy"])
print("Wet surface area:", hydrostatics["wet_surface_area"])
print("Displaced mass:", hydrostatics["disp_mass"])
print("Waterplane center:", hydrostatics["waterplane_center"])
print("Waterplane area:", hydrostatics["waterplane_area"])
print("Metacentric parameters:",
    body.transversal_metacentric_radius,
    body.longitudinal_metacentric_radius,
    body.transversal_metacentric_height,
    body.longitudinal_metacentric_height)

print(hydrostatics["hydrostatic_stiffness"])
# <xarray.DataArray 'hydrostatic_stiffness' (influenced_dof: 6, radiating_dof: 6)>
# [...]
# Coordinates:
#   * influenced_dof  (influenced_dof) <U7 'Surge' 'Sway' ... 'Yaw'
#   * radiating_dof   (radiating_dof) <U7 'Surge' 'Sway' ... 'Yaw'

print(hydrostatics["inertia_matrix"])
# <xarray.DataArray 'inertia_matrix' (influenced_dof: 6, radiating_dof: 6)>
# [...]
# Coordinates:
#   * influenced_dof  (influenced_dof) <U7 'Surge' 'Sway' ... 'Yaw'
#   * radiating_dof   (radiating_dof) <U7 'Surge' 'Sway' ... 'Yaw'

with open('Hydrostatic.txt','w') as fid:
    fid.write('Inertia\n')
    for i in range(6):
        tmp = hydrostatics["inertia_matrix"][i]
        np.savetxt(fid,[tmp],fmt='%.8e')
    fid.write('Hydrostatic stiffness\n')
    for i in range(6):
        tmp = hydrostatics["hydrostatic_stiffness"][i]
        np.savetxt(fid,[tmp],fmt='%.8e')

#Set up and solve problem
omega =0.501
wave_direction = np.pi/2
wave_amplitude= 1.0
radiation_problems = [cpt.RadiationProblem(body=body, radiating_dof=dof, omega=omega, water_depth=17.0, g=9.81, rho=1025) for dof in body.dofs]

diffraction_problem = cpt.DiffractionProblem(body=body, wave_direction=wave_direction, omega=omega,  water_depth=17.0, g=9.81, rho=1025)

solver = cpt.BEMSolver()
radiation_results = solver.solve_all(radiation_problems)

diffraction_result = solver.solve(diffraction_problem)
print(diffraction_result.forces)

dataset = cpt.assemble_dataset([diffraction_result] + radiation_results)
print(dataset["added_mass"])
print(dataset["radiation_damping"])
print(dataset["diffraction_force"])
print(dataset["Froude_Krylov_force"])
print(dataset["excitation_force"])

# # Calculations with a test matrix
theta_range = np.linspace(-np.pi/2,np.pi/2, 30)

# print(theta_range)
# Note: computation of Kochin function is activated with introduction of parameter 'theta'
# Nevertheless this computation fails if there is introduction of a lid mesh!

test_matrix = xr.Dataset(coords={
    'omega': [0.043, 0.044, 0.045, 0.046,   0.048,   0.049,   0.051,   0.052,   0.054,   0.056,
    0.058,   0.060,   0.062,   0.065,   0.068,   0.070,   0.074,   0.077,   0.081,   0.085,
    0.090,   0.095,   0.101,   0.108,   0.115,   0.124,   0.135,   0.147,   0.161,   0.179,
    0.200,   0.228,   0.265,   0.315,   0.388,   0.501,   0.697,   0.754,   0.819,   0.895,
    0.938,   0.985,   1.037,   1.095,   1.160,   1.235,   1.278,   1.324,   1.375,   1.432,
    1.496,   1.570,   1.655,   1.755,   1.877,   2.027,   2.221,   2.483,   2.867,   3.511],
    'wave_direction': [0.0,np.pi/6,np.pi/3,np.pi/2],
    # 'wave_direction': theta_range,
    'radiating_dof': list(body.dofs),
    'water_depth': [17.0],
    'rho': [1025],
    'theta': theta_range    
})


# test_matrix = xr.Dataset(coords={
#     'omega': [0.501,1.037],
#     'wave_direction': [np.pi/2],
#     'radiating_dof': list(body.dofs),
#     'water_depth': [17.0],
#     'rho': [1025],
#     'theta': [0.0,1.0]
# })
dataset = cpt.BEMSolver().fill_dataset(test_matrix, body)
rao = cpt.post_pro.rao(dataset)

# print(rao)
# print(dataset["added_mass"])
# print(dataset["radiation_damping"])
# print(dataset["diffraction_force"])
# print(dataset["Froude_Krylov_force"])
# print(dataset["excitation_force"])

num_freq = len(dataset["omega"])
num_wavedir = len(dataset["wave_direction"])
num_theta = len(dataset["theta"])
num_dof = len(dataset["radiating_dof"])

# Write added masses and radiation dampings in txt files
with open('Hydrodynamic.txt','w') as fid:
    for i in range(num_freq):
        fid.write('Wave frequency\n')
        tmp = dataset["omega"][i]
        np.savetxt(fid,[tmp],fmt='%.5e')
        fid.write('Added masses\n')
        for j in range(6):
            tmp = dataset["added_mass"][i][j][:]
            np.savetxt(fid,[tmp],fmt='%.8e')
    for i in range(num_freq):
        fid.write('Wave frequency\n')
        tmp = dataset["omega"][i]
        np.savetxt(fid,[tmp],fmt='%.5e')
        fid.write('Radiation damping\n')
        for j in range(6):
            tmp = dataset["radiation_damping"][i][j][:]
            np.savetxt(fid,[tmp],fmt='%.8e')

# Write the excitation forces and DAOs in txt files
with open('Forces.txt','w') as fid:
    for i in range(num_wavedir):
        fid.write('Wave direction\n')
        tmp = dataset["wave_direction"][i]
        np.savetxt(fid,[tmp],fmt='%.5e')
        fid.write(15*' '+'Frequency'+20*' '+'Surge'+30*' '+'Sway'+30*' '+'Heave'+30*' '+'Roll'+30*' '+'Pitch'+30*' '+'Yaw\n')
        for j in range(num_freq):
            tmp = dataset["omega"][j]
            tmp = np.append(tmp, dataset["excitation_force"][j][i][:])
            np.savetxt(fid,[tmp],fmt='%.8e')

with open('RAO.txt','w') as fid:
    for i in range(num_wavedir):
        fid.write('Wave direction\n')
        tmp = rao["wave_direction"][i]
        np.savetxt(fid,[tmp],fmt='%.5e')
        fid.write(4*' '+'Frequency'+6*' '+'Surge'+10*' '+'Sway'+11*' '+'Heave'+10*' '+'Roll'+11*' '+'Pitch'+10*' '+'Yaw\n')
        for j in range(num_freq):
            tmp = rao["omega"][j]
            tmp = np.append(tmp,abs(rao[j][i][:]))
            np.savetxt(fid,[tmp],fmt='%.8e')

# print(dataset['kochin'])
# print(dataset['kochin_diffraction'])

# Calculation of total Kochin function
wave_amplitude = 1   
# Diffraction component of the complex kochin function  
complex_kochin_diffraction = wave_amplitude * dataset['kochin_diffraction'] * np.exp(1j * (np.pi / 2))
# print(complex_kochin_diffraction)
# Radiation component of the complex kochin function
complex_kochin_radiation = np.zeros((num_freq,num_wavedir,num_theta)).astype(complex)
for i in range(num_freq):  
    for j in range(num_wavedir):  
        scalar = wave_amplitude *np.exp(-1j * (np.pi / 2)) * 1j * rao["omega"][i]
        array = np.zeros((num_freq,num_wavedir,num_theta)).astype(complex)
        # ind_theta used for determining the appropriate angle of kochin(w,dir,theta)
        ind_theta = -1
        for k in range(num_theta):
            if abs(dataset["theta"][k]-dataset["wave_direction"][j]) == np.min(abs(dataset["theta"]-dataset["wave_direction"][j])):
                ind_theta = k
                break
        if ind_theta != -1:
            print("element index  : ", ind_theta)
        else:
            print("The element not present")

        for l in range (num_dof):            
            array[i][j][ind_theta] = array[i][j][ind_theta] + dataset['kochin'][i][l][ind_theta]*rao[i][j][l]  
        complex_kochin_radiation[i][j][ind_theta] = scalar*array[i][j][ind_theta]

        # for l in range(num_theta):
        #     for k in range (num_dof):            
        #         array[i][j][l] = array[i][j][l] + dataset['kochin'][i][k][l]*rao[i][j][k]  
        #     complex_kochin_radiation[i][j][l] = scalar*array[i][j][l]

complex_kochin = np.add(complex_kochin_diffraction, complex_kochin_radiation)        

# print(complex_kochin_radiation)
# print(complex_kochin)

Fx,Fy = DF(dataset,complex_kochin,wave_amplitude,num_freq,num_wavedir,num_theta)

with open('Drift_Forces.txt','w') as fid:
    for i in range(num_wavedir):
        fid.write('Wave direction\n')
        tmp = dataset["wave_direction"][i]
        np.savetxt(fid,[tmp],fmt='%.5e')
        fid.write(1*' '+'Frequency'+10*' '+'Fx'+16*' '+'Fy\n')
        for j in range(num_freq):
            tmp = dataset["omega"][j]
            tmp = np.append(tmp, Fx[j][i])
            tmp = np.append(tmp, Fy[j][i])
            np.savetxt(fid,[tmp],fmt='%.8e')