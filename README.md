# Sim.B.A.D
A shipmooring program

The Sim.B.A.D. (SIMulation de Bateaux Amarrés sur Défenses) program allows the calculation of the mooring forces for a given type of vessel. PDSTRIP (Public Domain Strip Method) program, developed by H. Söding and V. Bertram provides the hydrodynamic parameters necessary to solve the linear hydrodynamic model in time domain (Cummins equation). It was simply adapted to obtain the characteristics of the desired parameters. 
This chain of programs is not intended to replace models that are more elaborate (3D models for example) but because of its quite fast operation, can provide useful information at the time of a pre-feasibility or feasibility study. The comparison with a certain number of results published in the literature (physical or mathematical models) shows that it can provide good orders of magnitude.
PDSTRIP:
PDSTRIP computes from the ship’s characteristics (cross-sectional strips of the wetted hull, displacement, position of centre of gravity, radii of gyration…) the following responses:
* Hydrostatic values, 
* Added masses,
* Radiation damping values,
* First order wave forces,
* Drift forces.
Pre-processing:
The hydrodynamic data pre-processing enables to provide:
* The “infinity” added masses, 
* The impulse response functions (“IRF”) of radiation damping values,
* The impulse response functions (“IRF”) of wave excitation and drift forces,
* Random wave time series according to different options:
o With random phase, on the basis of Hs, Tp and different spectrum density models (Jonswap, Bretschneider, Torsethaugen, McCormick, OchiHubble, Wallop) ,
o With the Fast method, on the basis of Hs, Tp and different spectrum density models (Jonswap, Bretschneider, Torsethaugen, McCormick, OchiHubble, Wallop) ,
o With random phase, on the basis of Hs, Tp and a specific density spectrum,
o Other options for wave time series: input as a specific file, regular waves.
* Time series of wave excitation forces,
* Time series of drift forces (3 options: no drift forces, mean forces, “Molin” variant).
SIMBAD program:
The SIMBAD program can be used to simulate the behavior of a moored ship, using the Cummins equation (time simulation) and taking into account the following elements:
* Random wave forces (excitation and drift, see above),
* Gusty wind forces (generation of wind speeds according to a wind spectrum),
* Forces induced by a passing ship,
* Linear or non-linear mooring characteristics for mooring lines and fenders, as well as rigid connection type mooring system.
Post-processing:
The results data post-processing enables to provide plots of time series and analysis of results (min, max, 1/3, 1/10 values) for:
* Mooring line and fender loads,
* Motions, speeds and accelerations at the center of gravity and at other possible points of the ship.

