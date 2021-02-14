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

Milestones: Positions, velocities and energy are correctly stored in their respective arrays. The functions for the forces, Euler method and energy have been implemented, but we are still missing the periodic boundary conditions, and need to figure out the forces so that they can actually make an impact on the system.

(due 14 February 2021, 23:59)


## Week 2
(due 21 February 2021, 23:59)


## Week 3
(due 28 February 2021, 23:59)


## Week 4
(due 7 March 2021, 23:59)


## Week 5
(due 14 March 2021, 23:59)
