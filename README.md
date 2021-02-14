# Project 1: Molecular dynamics simulation of Argon atoms

The README.md file serves as a reference for other users visiting your repository.
It should contain a brief description of your project, and document the steps others need to take to get your application up and running.
In addition, it should list the authors of the project.

Authors:               Student Numbers
Alberto Gori           (5391776)
Matteo de Luca         (5388783)
Pol de Dalmau Huguet   (5414024)


How to use our code:
For now only look at the Week1_complete.ipynb notebook. Here we have written a script that canimplement the Euler Method almost as intended. One can adjust the number of time steps, time steps, and mass of the argon atoms. Other variables are defined but changing them won't have the intended effect. The output will be at the bottom (cell with for-loop) consisting of two matrices:
final_matrix_pos: [[x1, y1, z1],
                   [x2, y2, z2],
                   [x3, y3, z3]]
final_matrix_vel: [[v_x1, v_y1, v_z1],
                   [v_x2, v_y2, v_z2],
                   [v_x3, v_y3, v_z3]]

Also the total energy is calculated at every time step. No visualization is yet available. Also, the boundary condition is only applied to the position. The nearest distance to calculate forces is not yet implemented. We only use the distance within the box.