# Weekly progress journal

## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on the day of the lecture** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented.
2. write about your progress **until Sunday, 23:59** before the next lecture with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj1-moldyn-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly
useful.

## Week 1
10.02.2021: Decided an effective method to write and test our code is to use the notebook directly. First, we agreed on what types of functions we'll make and which inputs and outputs they have so we can make them compatible with one another. 

We decided the following distribution of work:
- Matteo: Function that calculates the forces on the particles (input: positions of all particles at all times; output: Force components at most recent timestep)
- Alberto: Finish the function that stores each particle's velocity and position at every timestep (input: init_pos and init_vel; output: longer matrix with all coordinates at all timesteps)
and implement boundary condition
- Pol: Make a function that implements Euler's Method (input: velocities and positions up to timestep n, output: positions and velocities at most recent timestep n+1)

Milestones: 

- Positions, velocities and energy are correctly stored in their respective arrays [Skeleton lines 323-330](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L323-330).
- The functions for the forces [(Forces: Skeleton lines 165-196)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L165-196), Euler method [(Euler: Skeleton lines 200-206)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L200-206), kinetic energy [(Kinetic: Skeleton lines 219-234)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L219-234) and potential energy [(Potential: Skeleton lines 237-259)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L237-259) have been implemented, but we are still missing the periodic boundary conditions, and need to figure out the forces so that they can actually make an impact on the system.

(due 14 February 2021, 23:59)


## Week 2
17.02.21: We worked together on the dimensionless units and the kinetic energy expression, and then divided the work: 
- Matteo: Improve the forces function to get rid of for loops, as well as implement it for a 3D system
- Alberto: Improve the journal and plot the distances and energies
- Pol: Create the atomic_distances function and switch the code to 3D

Milestones: 

- We were able to implement the dimensionless units and compute the energy using the minimal image convention implemented in the atomic_distances function [(Atomic distances: Skeleton lines 126-158)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L126-158) and we were able to plot the relative distances and the energies over time. After a few steps, the relative distance reaches a minimum (close to zero) for which the kinetic energy becomes massive and the total energy is not conserved.
Compared to week 1, we removed the majority of useless for loops and replaced them with numpy funcionalities, but some for loops are still there. 
- LJ forces without loops: [(LJ: Skeleton lines 165-196)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L165-196)
- Energy without loops: kinetic [(Kinetic: Skeleton lines 219-234)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L219-234) and potential [(Potential: Skeleton lines 237-259)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L237-259)
(due 21 February 2021, 23:59)


## Week 3
- Managed to correct errors with previous week's code (see Week2_complete.ipynb) and implement verlet method [(Verlet: Skeleton lines 208-215)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L208-215) for more than 2 particles (only tested for 3 particles thus far). Many of the for loops that slow down the code have not been removed. No attempt at structuring the code has been done this week.
- The functions calculating the potential energy [(Potential: Skeleton lines 237-259)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L237-259) and the lj force [(LJ: Skeleton lines 165-196)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L165-196) were promoted to n-particles.


(due 28 February 2021, 23:59)


## Week 4
Please refer to the Week4 velocities and rescaling.ipynb (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week4%20velocities%20and%20rescaling.ipynb#L361) to see a successful implementation of the rescaling. The fcc lattice and appropriate velocity initializations are still in progress.

Edit: this week we were not able to work on the code properly due to other obligations.
(due 7 March 2021, 23:59)


## Week 5
- Pol: implemented the fcc_latice function.
- Alberto and Matteo: worked on the Maxwell-Boltmzann distribution of velocities and rescaling.
- Matteo: implemented error function.
- Alberto: implemented pressure function and restructured the code.
- The code works, is now structured in a better way and each function has its description, but the results are not always reasonable. 
- We successfully implemented the Maxwell-Boltzmann distribution of velocities [(M-B: Skeleton lines 15-16)](https://gitlab.kwant-project.org/computational_physics/projects/Project-  1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L15-16) and plotted it [Skeleton lines 24-27](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L24-27) (Figure 4 in figures folder). 
- The fcc lattice function works and we successfully plotted it in 3-D (Figure 5 in figures folder). [(FCC: Skeleton lines 29-123)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L29-123)
- The pressure function works, but the resulting pressures are sometimes huge. [(Pressure: Skeleton lines 262-273)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L262-273) gives the value that will be averaged, while in the simulate function [(Pressure: Skeleton lines 342-345)](https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/skeleton.py#L342-345) we compute the pressure averaging.

The energies plot for 16 particles at 10K with 10000 timesteps of 0.001 is shown in figure 7. It is probably our best result so far.
Figure 8 also shows the error and the energies for 100K while the other parametrs are the same as in figure 7. The pressure was 4.97 in dimensionless units.

To valide our simulation, we will compute the pressure and compare it to literature values. Also, simulating the code with certain initial values for the temparature could give us information on eventual changes of state. 

With the actual performance of our code we can run simulations with number of timesteps of the order of $10^4$ and a variable number of particles. The number of particles is the major factor in computing time. With a number of particles between 10 and 20 (and 10000 timesteps) it takes a few minutes for the entire simulation to run.


Questions:
- The total energy keeps oscillating along with the potential one, even for 2 particles. In some configurations, the potential and total energy keep getting smaller while the kinetic energy is almost zero. The potential energy is very negative, up to -140 for 16 particles in a 3x3x3 box. The initial potential energy in this case was already -128.  Even for two particles with zero initial velocity and without rescaling, the total energy oscillates with the potential (please see figure 1 in the folder for plot).
- Pressure is huge even for big boxes compared to the number of particles. For 16 particles with a 3x3x3 box at 10K (0.0001 time step duration, 10000 time steps) the pressure is 63.25 in dimensionless units, which is around 26000 bars, definitely too big (plot is in figure 2 inside the folder). To switch from dimensionless units to SI units for pressure, we multiplied by $(\sigma^3)/\epsilon$ .
- Rescaling works, but we are not sure that it is correct (please see figure folder for plots).
(due 14 March 2021, 23:59)
