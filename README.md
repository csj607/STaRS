# STaRS
Sejong Radiative Transfer through Raman and Rayleigh Scattering with atomic hydrogen

In Astronphysics, Raman spectroscopy is good tool to investigate Symbiotic Stars, Planetary Nebulae, and Active Galacitive Nuclei.

STaRS is the code for Radiative Transfer through Raman and Rayleigh Scattering with atomic hydrogen.
This code is 3D grid based Monte Carlo simulation tracing each generating photon packet.
The information of the photon packet include wavelength, position, and polarization.

The basic langauge and compiler are FORTRAN and intel FORTRAN.
I adopted parallel computing and shared memory technique for fast calculating and handling a memory.
If you have any question about the code, you send the email to "csj607@gmail.com".
Any comments for development and suggestions for collaboration are well come.

The paper for STaRS is accepted for publication in JKAS.
The link for the manuscript in astro-ph,
https://arxiv.org/abs/2012.03424

The title is "3D Grid-Based Monte Carlo Code for Radiative Transfer through Raman and Rayleigh Scattering
with Atomic Hydrogen --- STaRS".

Publication List using STaRS
Bo-Eun Choi, Hee-Won Lee, 2020, ApJL, 903, L39, Discovery of Raman-scattered He II λ6545 in the Planetary Nebulae NGC 6886 and NGC 6881
Bo-Eun Choi, Seok-Jun Chang, Ho-Gyu Lee, and Hee-Won Lee, 2020, ApJ, 899, 12C, Line Formation of Raman-scattered He II λ 4851 in an Expanding Spherical H I Shell in Young Planetary Nebulae

Source files

main.f90 : the main code to run STaRS

RT_grid.f90 : the module to set the scattering geometry

RT_photon.f90 : the module to generate the initial photons and describe the scattering process of photon. 

RT_obs.f90 : the module to collect the information of escaping photon.

RT_cross.f90 : the module to compute the scattering cross section and braching ratio by atomic physic.

random_mt.f90 : the module including the random generator

memory_mod.f90 : the module for the shared memory technique.

HOW TO USE THE CODE

1. Download all of file .f90 (source code) and com.sh (commands)
2. chmod +x com.sh
3. Run com.sh
3. Run STaRS.run file using the command, 'mpirun'
