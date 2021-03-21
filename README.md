# Project 1: Molecular dynamics simulation of Argon atoms

The README.md file serves as a reference for other users visiting your repository.
It should contain a brief description of your project, and document the steps others need to take to get your application up and running.
In addition, it should list the authors of the project.

Authors:               Student Numbers
Alberto Gori           (5391776)
Matteo de Luca         (5388783)
Pol de Dalmau Huguet   (5414024)


How to use our code:
Open the Test.ipynb notebook and run each cell from top to bottom. Here we have written a script that can implement the Verlet Method in 3D for more n particles. One can adjust the number of time steps, number of particles, time step size, system size L, the lattice constant a of the fcc lattice and the temperature. One will obtain a plot of energy (kinetic, potential and total) which is rescaled several times and a value for the pressure in normalized units. Also, if the boolean values in the simulate function are set to "True", the system will also return a pair correlation function and the error on the pressure.
