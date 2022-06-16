# Project 1: Molecular dynamics simulation of Argon atoms

Authors:               Student Numbers
Alberto Gori           (5391776)
Matteo de Luca         (5388783)
Pol de Dalmau Huguet   (5414024)


How to use our code:
Open the Main.ipynb notebook and run each cell from top to bottom. Here we have written a script that can implement the Verlet Method in 3D for n particles. One can adjust the number of time steps, number of particles, time step size, system size L (volume is LxLxL), of the fcc lattice and the temperature in natural units. The density of particles must be inferred from the choice of n and L. After simulating, the outputs are vectors for kinetic, potential and total enery, pair correlation values, pressure, error of the pressure using data blocking and the density. Then the functions to plot the energies, error and pair correlation function can be called after running a simulation to visualize the results.
