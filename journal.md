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
Matteo: Function that calculates the forces on the particles (input: positions of all particles at all times; output: Force components at most recent timestep)
Alberto: Finish the function that stores each particle's velocity and position at every timestep (input: init_pos and init_vel; output: longer matrix with all coordinates at all timesteps)
and implement boundary condition
Pol: Make a function that implements Euler's Method (input: velocities and positions up to timestep n, output: positions and velocities at most recent timestep n+1)

Milestones: Positions, velocities and energy are correctly stored in their respective arrays (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week1_complete.ipynb#L385-387).
The functions for the forces (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week1_complete.ipynb#L60-84), Euler method (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week1_complete.ipynb#L97-104) and energy (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week1_complete.ipynb#L114-142) have been implemented, but we are still missing the periodic boundary conditions, and need to figure out the forces so that they can actually make an impact on the system.

(due 14 February 2021, 23:59)


## Week 2
17.02.21: We worked together on the dimensionless units and the kinetic energy expression, and then divided the work: 
Matteo: Improve the forces function to get rid of for loops, as well as implement it for a 3D system
Alberto: Improve the journal and plot the distances and energies
Pol: Create the atomic_distances function and switch the code to 3D

Milestones: we were able to implement the dimensionless units and compute the energy using the minimal image convention implemented in the atomic_distances function (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week2_complete.ipynb#L157-178) and we were able to plot the relative distances and the energies over time (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Plot%20of%20positions%20and%20energies%20-%20Week2.png). After a few steps, the relative distance reaches a minimum (close to zero) for which the kinetic energy becomes massive and the total energy is not conserved.
Compared to week 1, we removed the majority of useless for loops and replaced them with numpy funcionalities, but some for loops are still there. 
LJ forces without loops: https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Plot%20of%20positions%20and%20energies%20-%20Week2.png
Energy without loops: kinetic (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week2_complete.ipynb#L228-231) and potential (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week2_complete.ipynb#L228-231)
(due 21 February 2021, 23:59)


## Week 3
Managed to correct errors with previous week's code (see Week2_complete.ipynb) and implement verlet method (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week3.ipynb#L187-194) for more than 2 particles (only tested for 3 particles thus far). Many of the for loops that slow down the code have not been removed. No attempt at structuring the code has been done this week.
The functions calculating the potential energy (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week3.ipynb#L229-239) and the force (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week3.ipynb#L146-165) were promoted to n-particles.


(due 28 February 2021, 23:59)


## Week 4
Please refer to the Week4 velocities and rescaling.ipynb (https://gitlab.kwant-project.org/computational_physics/projects/Project-1_albertogori_compphys_bot_matteodeluca_pdedalmauhugue/-/blob/master/Week4%20velocities%20and%20rescaling.ipynb#L361) to see a successful implementation of the rescaling. The fcc lattice and appropriate velocity initializations are still in progress.
(due 7 March 2021, 23:59)


## Week 5
(due 14 March 2021, 23:59)
