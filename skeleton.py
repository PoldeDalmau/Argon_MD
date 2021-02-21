"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import numpy as np
import matplotlib.pyplot as plt

n = 2   #number of particles
L = 1   #System is of size LxL with the origin at the center. The units of L are of 1/s
init_pos = np.array([[0.1, 0.8, 0.8],[0.7, 0.9, 0.9]])
init_vel = np.random.uniform(low=-1, high=1, size=(n, 3))
e =  1.65 * pow(10,-21)
s =  3.405 * pow(10,-10)
time_step = 1e-3
m =  6.6335 * pow(10, -26) #mass of argon atom in kg
number_of_steps = 10
dimensions = 3
times = np.arange(0,(number_of_steps + 1) * time_step, time_step)

final_matrix_pos = np.copy(init_pos)
final_matrix_vel = np.copy(init_vel)

#the function below computes the atomic distances between particles
W = 1   #Width 
K = 1   #Depth
def atomic_distances(pos, output, L=1, W=1): # output = 0 gives the relative distances, output = 1 gives the relative positions
    rel_pos = pos[:, None, :] - pos[None, :, :]                 # returns one matrix for each particle. Relative distances within the box
    rel_dist = np.zeros((n, n, 1))
    for i in range(n):
        for k in range(n):
            for l in range(dimensions):
                if l == 0: # x-values i.e. change by L
                    rel_pos[i,k,l] = min(rel_pos[i,k,l], (rel_pos[i,k,l]+L)%L, (rel_pos[i,k,l]-L)%L, rel_pos[i,k,l]+L, rel_pos[i,k,l]-L, key=abs) # takes the smallest distance comparing all images
                if l == 1:
                    rel_pos[i,k,l] = min(rel_pos[i,k,l], (rel_pos[i,k,l]+W)%W, (rel_pos[i,k,l]-W)%W, rel_pos[i,k,l]+W, rel_pos[i,k,l]-W, key=abs)
                if l == 2:
                    rel_pos[i,k,l] = min(rel_pos[i,k,l], (rel_pos[i,k,l]+K)%K, (rel_pos[i,k,l]-K)%K, rel_pos[i,k,l]+K, rel_pos[i,k,l]-K, key=abs)
    for i in range(n):    
        for k in range(n):
            rel_dist[i,k] = np.sqrt(np.dot(rel_pos[i,k,:],rel_pos[i,k,:]))
    
    if output == 0:
        return rel_dist
    elif output == 1:
        return rel_pos
    else:
        print("The second parameter in the lj_forces must be either 0(for distances) or 1(for positions)")
        
#compute the forces on the particles at each timestep, pos = [m]/[s] = 
def lj_force(position, n):
    
    rel_pos = atomic_distances(position[:,-3:], 1) # with pos = np.array(N_particels, 3)
    r = np.sqrt(np.sum((rel_pos)**2, axis=2)) #n x n simmetric matrix, r[i,j] is the distance between the i-th and the j-th particles
    r = r[0,1]
    F = np.zeros((n,3))
    F = 24*rel_pos*(((1/r)**12)-((1/r)**6))
    F_matrix=np.sum(F,axis=0)    
    
    return F_matrix

#the function below implements the euler method
def euler(final_matrix_pos, final_matrix_vel):
    
    latest_pos = np.copy(final_matrix_pos[:,-3:])         # -3 is for the dimensions
    latest_vel = np.copy(final_matrix_vel[:,-3:])
    
    new_latest_pos = latest_pos + latest_vel * time_step
    new_latest_vel = latest_vel + lj_force(latest_pos, n)
    return new_latest_pos, new_latest_vel

#function to compute the kinetic energy
def kin_en(v): #v is the last step velocity
    K = 0.5*np.sum(v**2)
    return K

#function to compute the potential energy
def pot_en(position, n): #position is the matrix with all the positions stored in it, n is the number of particles
    current_pos = np.copy(position[:,-3:])
    rel_pos = atomic_distances(current_pos, 1) # with pos = np.array(N_particels, 3)
    r = np.sqrt(np.sum((rel_pos)**2, axis=2)) #n x n simmetric matrix, r[i,j] is the distance between the i-th and the j-th particles
    r = r[0,1]
    U = 4*(1/r**12 - 1/r**6)
    
    return U

def simulate():
    #Create a 2x8 matrix to store the velocity of each particle at each step in time.
    next_step_velocity = np.copy(init_vel)

    #Create a 2x8 matrix to store the position of each particle at each step in time.
    next_step_position = np.copy(init_pos)

    final_vector_kin = np.array([kin_en(init_vel)])
    final_vector_pot = np.array([pot_en(final_matrix_pos, n)])
    final_vector_energy = np.array([kin_en(init_vel) + pot_en(final_matrix_pos, n)])

    final_rel_dist = atomic_distances(init_pos, 0)

    for i in range(number_of_steps):
        next_step_position, next_step_velocity = euler(final_matrix_pos, final_matrix_vel)
        for k in range(n): #implement boundary conditions
            for j in range(n):
                next_step_position[k, j] = next_step_position[k, j]% L
        final_matrix_pos =  np.concatenate((final_matrix_pos, next_step_position), axis=1, out=None)
        final_matrix_vel =  np.concatenate((final_matrix_vel, next_step_velocity), axis=1, out=None)
        final_vector_kin = np.concatenate((final_vector_kin, np.array([kin_en(next_step_velocity)])), axis = 0, out = None)
        final_vector_pot = np.concatenate((final_vector_pot, np.array([pot_en(final_matrix_pos, n)])), axis=0, out=None)
        final_vector_energy = np.concatenate((final_vector_energy, np.array([kin_en(next_step_velocity) + pot_en(final_matrix_pos, n)])), axis=0, out=None)
        final_rel_dist = np.concatenate((final_rel_dist, atomic_distances(next_step_position, 0)), axis = 1, out = None)

    #print("Positions:\n" , final_matrix_pos)
    #print("Velocities:\n" , final_matrix_vel)
    #print("Energy:\n" , final_vector_energy)
    #print(final_vector_kin)
    #print(final_vector_pot)
    #print(final_rel_dist)


    plt.plot(times,final_rel_dist[0, 1::2, ])
    plt.title("Relative distances over time")
    plt.xlabel("$t/(m\sigma^2/\epsilon)^{1/2}$")
    plt.ylabel("$r/\sigma$")
    plt.show()

    plt.plot(times, final_vector_energy)
    plt.plot(times, final_vector_kin)
    plt.plot(times, final_vector_pot)
    plt.legend(["Total Energy", "Kinetic Energy", "Potential Energy"])
    plt.title("Energy over time")
    plt.xlabel("$t/(m\sigma^2/\epsilon)^{1/2}$")
    plt.ylabel("$E/\epsilon$")
    plt.show()

#we are keeping the initial functions to have the description available at all times
def simulate_empty(init_pos, init_vel, num_tsteps, timestep, box_dim): #arguments are global

    """
    Molecular dynamics simulation using the Euler or Verlet's algorithms
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep. HELOOOOOO

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    Any quantities or observables that you wish to study.
    """
    
    
    
    return 


def atomic_distances_empty(pos, box_dim):
    """
    Calculates relative positions and distances between particles.

    parameters
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_dim : float
        The dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles
    rel_dist : np.ndarray
        The distance between particles
    """

    return


def lj_force_empty(rel_pos, rel_dist):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    np.ndarray
        The net force acting on particle i due to all other particles
    """

    return


def fcc_lattice(num_atoms, lat_const):
    
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system
    lattice_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """

    return


def kinetic_energy_empty(vel):
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    float
        The total kinetic energy of the system.
    """

    return


def potential_energy_empty(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    float
        The total potential energy of the system.
    """

    return


def init_velocity_empty(num_atoms, temp):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities
    """

    return

#our code starts here!

import numpy as np
import matplotlib.pyplot as plt



W = 1   #Width 
K = 1   #Depth
def atomic_distances(pos, output, L=1, W=1): # output = 0 gives the relative distances, output = 1 gives the relative positions
    rel_pos = pos[:, None, :] - pos[None, :, :]                 # returns one matrix for each particle. Relative distances within the box
    rel_dist = np.zeros((n, n, 1))
    for i in range(n):
        for k in range(n):
            for l in range(dimensions):
                if l == 0: # x-values i.e. change by L
                    rel_pos[i,k,l] = min(rel_pos[i,k,l], (rel_pos[i,k,l]+L)%L, (rel_pos[i,k,l]-L)%L, rel_pos[i,k,l]+L, rel_pos[i,k,l]-L, key=abs) # takes the smallest distance comparing all images
                if l == 1:
                    rel_pos[i,k,l] = min(rel_pos[i,k,l], (rel_pos[i,k,l]+W)%W, (rel_pos[i,k,l]-W)%W, rel_pos[i,k,l]+W, rel_pos[i,k,l]-W, key=abs)
                if l == 2:
                    rel_pos[i,k,l] = min(rel_pos[i,k,l], (rel_pos[i,k,l]+K)%K, (rel_pos[i,k,l]-K)%K, rel_pos[i,k,l]+K, rel_pos[i,k,l]-K, key=abs)
    for i in range(n):    
        for k in range(n):
            rel_dist[i,k] = np.sqrt(np.dot(rel_pos[i,k,:],rel_pos[i,k,:]))
    
    if output == 0:
        return rel_dist
    elif output == 1:
        return rel_pos
    else:
        print("The second parameter in the lj_forces must be either 0(for distances) or 1(for positions)")