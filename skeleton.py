"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""


import numpy as np
import matplotlib.pyplot as plt


#n = 3   #number of particles
#L = 1000   #System is of size LxL with the origin at the center. The units of L are of 1/s
#init_pos = np.array([[0, 0, 0],[L-1, L-1, L-1], [1, 1, 1]])
#init_vel = np.random.normal(np.sqrt(T/119.8), size = (n, 3))
#init_vel = init_vel - np.mean(init_vel, axis=0)
#e =  1.65 * pow(10,-21)
#s =  3.405 * pow(10,-10)
#time_step = 1e-3
#m =  6.6335 * pow(10, -26) #mass of argon atom in kg
#number_of_steps = 10000
T = 1.5 #initial temperature of the system
times = np.arange(0,(number_of_steps + 1) * time_step, time_step)
#test for Boltzmann distribution. Must be run with at least 1000 particles
#v_norm = np.linalg.norm(init_vel, axis=1)
#plt.hist(v_norm, 20)
#plt.show()

def cleaner(init_pos, strict, box_size = L):
    """checks all values greater than the box size (outside of fcc unit cell) and deletes them
    
    Parameters
    ----------
    init_pos : np array (some number, 3)
        configuration of particles
    strict : 
        takes 2 possible values. determines whether it muust remove particles on and beyond the boundary ("very")
        or only the ones beyond the boundary ("less")
    Returns
    -------
    init_pos : np.ndarray
        Array of particle coordinates
    """
    guilt = []          #list to be filled with all column numbers of forbidden points"
    for i in range(len(init_pos)):
            for k in range(3):
                if (init_pos[i,k] >= box_size or init_pos[i,k] < 0) and strict == "very":
                    guilt.append(i)
                if (init_pos[i,k] > box_size or init_pos[i,k] < 0) and strict == "less":
                    guilt.append(i)
    #remove elements counted twice
    res = [] 
    [res.append(x) for x in guilt if x not in res]  
    init_pos = np.delete(init_pos, res[::-1],0)                      # very important to start with the last one and go backwards
    
    return (init_pos)
    #return "unit_sc", unit_sc.shape, "unit_sc_1",  unit_sc_1.shape, "init_pos", init_pos.shape
    #return "unit_sc", unit_sc, "unit_sc_1",  unit_sc_1, "init_pos", init_pos
    
def fcc_lattice(a):
    
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    a : float
        The lattice constant for an fcc lattice.

    Returns
    -------
    init_pos : np.ndarray
        Array of particle coordinates
    """
    square_xy = np.array([[0, 0, 0], [a, 0, 0], [a, a, 0],[0, a, 0]])
    alongx = np.copy(square_xy)
    # This will begin to put squares one next to another 
    # along the x-axis. Distance between squares must be 2a
    multiples = np.arange(0, L, 2*a)
    multiples_top = np.arange(0, L, a)
    
    multiples = np.reshape(multiples, (len(multiples),1))
    multiples_top = np.reshape(multiples_top, (len(multiples_top),1))
    
    multiples_x = multiples*np.array([1, 0, 0])
    multiples_y = multiples*np.array([0, 1, 0])
    multiples_z = multiples_top * np.array([0, 0, 1])
    
    for i in range(1, len(multiples)):
        to_add = multiples_x[i] + square_xy
        alongx = np.concatenate((alongx,to_add), axis = 0) 
    alongx = cleaner(alongx, "less")
    

    xy_plane = alongx  
    for i in range(1, len(multiples)):
        to_add1 = multiples_y[i] + alongx
        xy_plane = np.concatenate((xy_plane,to_add1), axis = 0)
    xy_plane = cleaner(xy_plane, "less")
    
    sc_lattice = xy_plane
    for i in range(1, len(multiples_top)):                                 
        to_add2 = multiples_z[i] + xy_plane
        sc_lattice = np.concatenate((sc_lattice,to_add2), axis = 0)
    sc_lattice = cleaner(sc_lattice, "less")

# an fcc is equivalent to four simple cubic lattices,
    sc_1 = np.copy(sc_lattice) + a/2*np.array([1,1,0])                         #   one with a(1/2,1/2,0) offset, 
    sc_2 = np.copy(sc_lattice) + a/2*np.array([0,1,1])                         #   one with a(0,1/2,1/2) offset,              
    sc_3 = np.copy(sc_lattice) + a/2*np.array([1,0,1])                         #   and one with a(1/2,0,1/2) offset           
    fcc_positions = np.concatenate((sc_lattice, sc_1, sc_2, sc_3), axis = 0)
    
    
    fcc_positions = np.concatenate((fcc_positions, fcc_positions + L/2 * np.array([0, 0, 1])), axis = 0)     # corrects the fact we were only making half the particles we were supposed to
    fcc_positions = cleaner(fcc_positions, "very")
    max_n = int((L/a)**3 * 4)            #Theoretical maximum of particles
    if n > max_n:
        print("You wanted %i particles, but only %i particles can actually fit here!" % (n, max_n))
        fcc_positions = fcc_positions[:max_n,:] #picks n particles out of all, 
   
    fcc_positions = fcc_positions[:n,:]  #picks n particles out of all,     
    return (fcc_positions)

def init_velocity(T, n):
    """
    Initializes the initial velocities of the particles.
    
    Parameters
    ----------
    T : float
        Initial temperature of the system
    n : integer
        Number of particles inside the system
    
    returns
    -------
    init_vel : np.ndarray
        nx3 array containing the initial velocities of the particles. Each column represents an axis, each row a different particle
    """
    init_vel = np.zeros((n,3))
    init_vel = np.random.normal(np.sqrt(T/119.8), size = (n, 3))
    #init_vel = init_vel - np.mean(init_vel, axis=0)
    
    return init_vel


def atomic_distances(pos, output): # output = 0 gives the relative distances, output = 1 gives the relative positions
    """
    Calculates relative positions and distances between particles.

    parameters
    pos : np.ndarray
        The positions of the particles in cartesian space at the latest timestep.
    output: float
        If output = 0 the function returns the relative distances between the particles.
        If output = 1 the function returns the relative positions of the particles. 
        If output is neither 1 or 0, an error message is printed.

    returns
    -------
    r : np.ndarray
        Relative distances between the particles.
    rel_pos : np.ndarray
        The relative positions of the particles.
    """
    rel_pos = pos[:, None, :] - pos[None, :, :]                 # returns one matrix for each particle. Relative distances Lithin the box
    n = np.shape(rel_pos)[0] # number of particles
    rel_pos = ((rel_pos + L/2) % L) - L/2              
    r = np.sqrt(np.sum((rel_pos)**2, axis=2)) # n x n simmetric matrix, r[i,j] is the distance between the i-th and the j-th particles
    
    if output == 0:
        return r
    elif output == 1:
        return rel_pos
    else:
        print("The second parameter in the atomic_distances must be either 0(for distances) or 1(for positions)")
        
        
        
        
#compute the forces on the particles at each timestep

def lj_force(rel_pos, r):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions at the latest timestep.
    r : np.ndarray
        Relative particle distances at the latest timestep.
    

    Returns
    -------
    np.ndarray
        The net force acting on particle i due to all other particles.
    """
    #rel_pos = atomic_distances(position[:,-3:], 1) # calculates the relative position of the particles, with pos = np.array(N_particels, 3)
    #r = atomic_distances(position[:,-3:], 0) # n x n simmetric matrix, r[i,j] is the distance between the i-th and the j-th particles
    zeros = np.zeros((n,n))
    np.fill_diagonal(zeros, np.inf)         # nxn diagonal matrix with np.inf on the diagonal
    r = r + zeros
    temp = np.reshape(r, (n,n,1))           # temporary array to build R in the right shape (n,n,3)
    R = np.concatenate((temp, temp, temp), axis=2)
    # R are n matricies in the form nx3, where the i-th matrix contains the distances from the i-th particle, calling R_i the i-th matrix, 
    # R_i[j,k] gives the distance between the i-th and the j-th particles and is the same value for all the 'k in range(3)', i.e. on all the j-th line
    # in addition, R_i[i,k] (i-th line of the i-th matrix) is equal to 1 instead of 0 in order not to divide by 0 in the computation of F
    # we can do this since rel_pos_i[i,k] multiplies everything in F and is equal to 0
    F = np.zeros((n,n,3))
    F = -24*rel_pos*((2*(1/R)**12)-((1/R)**6)) 
    F_matrix = np.sum(F, axis=0)  
    return F_matrix # The output is an nx3 matrix
    
    
    
def euler(latest_pos, latest_vel, rel_pos, rel_dist):
    """
    Updates the positions of each particle at each timestep using Euler algorithm.
    
    Parameters
    ----------
    latest_pos : np.ndarray
        Matrix with the positions of the particles at the latest timestep
    latest_vel : np.ndarray
        Matrix with the velocities of the particles at the latest timestep
        
    Return
    ------
    np.ndarray
        Returns two arrays: one containing the postions, and one containing the velocities at the latest timestep
    """
    #latest_pos = np.copy(final_matrix_pos[:,-3:])# -3 is for the dimensions
    #latest_vel = np.copy(final_matrix_vel[:,-3:])
    
    new_latest_pos = latest_pos + latest_vel * time_step
    new_latest_vel = latest_vel + time_step * lj_force(rel_pos, rel_dist)
    return new_latest_pos, new_latest_vel

def verlet(latest_pos, latest_vel, rel_pos, rel_dist):
    """
    Updates the positions of each particle at each timestep using Verlet algorithm.
    
    Parameters
    ----------
    latest_pos : np.ndarray
        Matrix with the positions of the particles at the latest timestep
    latest_vel : np.ndarray
        Matrix with the velocities of the particles at the latest timestep
        
    Return
    ------
    np.ndarray
        Returns two arrays: one containing the postions, and one containing the velocities at the latest timestep
    """
    
    #latest_pos = np.copy(final_matrix_pos[:,-3:])# -3 is for the dimensions
    #latest_vel = np.copy(final_matrix_vel[:,-3:])
    
    new_latest_pos = latest_pos + time_step * latest_vel + (1/2) * time_step**2  * lj_force(rel_pos, rel_dist)
    new_latest_vel = latest_vel + (time_step/2) * (lj_force(rel_pos, rel_dist) + lj_force(rel_pos, rel_dist))
    return new_latest_pos, new_latest_vel


#function to compute the kinetic energy
def kin_en(v): #v is the last step velocity
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    v: np.ndarray
        Velocity of particles at the latest timestep.

    Returns
    -------
    float
        The total kinetic energy of the system.
    """
    K = 0.5*np.sum(v**2)
    return K
    
#function to compute the potential energy
def pot_en(r): #position is the last step position of the particles
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    r : np.ndarray
        Particles relative distances at the latest timestep.

    Returns
    -------
    float
        The total potential energy of the system.
    """
    #rel_pos = atomic_distances(position, 1) # with pos = np.array(n, 3)
    #r = atomic_distances(position, 0) #n x n simmetric matrix, r[i,j] is the distance between the i-th and the j-th particles
    i = np.where(r == 0) # indices where the r matrix is 0
    ones = np.zeros((n,n)) # nxn matrix
    ones[i] = 1 # matrix with ones where r has zeros
    R = r + ones  # same matrix as r, but with ones where r has zeros. we do this to avoid dividing by zero
    u = 1/2*4*(1/R**12 - 1/R**6) # factor 1/2 to compensate for double counting
    U = np.sum(u)
    return U

#function to compute the pressure
def pressure(r): 
    """
    Computes the value that will be averaged to find the pressure
    
    Parameters
    ----------
    r : np.ndarray
        Particles relative distances at the latest timestep.

    Returns
    -------
    float
        Returns the sum of the element of the matrix, which is the value that will be averaged to find the pressure.
    """
    
    #r = atomic_distances(position[:,-3:], 0) #n x n simmetric matrix, r[i,j] is the distance between the i-th and the j-th particles
    i = np.where(r == 0) # indices where the r matrix is 0
    inf = np.zeros((n,n)) # nxn matrix
    inf[i] = np.inf # matrix with inifinity where r has zeros
    R = r + inf  # same matrix as r, but with ones where r has zeros. we do this to avoid dividing by zero
    tba_matrix = -12*((2*(1/R)**12)-((1/R)**6))
    tba = np.sum(tba_matrix)
    #rho = n/(L*L*L)
    #P = (rho*119.8/T) - (rho/(3*n))*
    return tba

def pair_correlation(rel_dist, dr):
    """
    Calculates pair correlation for all particles at a given time.
    Parameters
    ----------
    rel_dist : np.ndarray
        relative distances between atoms as calculated by atomic_distances()
    Returns
    -------
    Fills a histogram.
        
    """

    rel_dist = np.reshape(rel_dist, (n*n))
    bins_list = np.arange(dr, L+dr, dr)
    counts = pd.cut(rel_dist, bins_list, include_lowest = True)
    counts= counts.value_counts()
    counts = counts.to_numpy()
    return counts/2 #avoids double counting
    
    
def simulate(algorithm, rescaling_bool, pressure_bool, error_bool, pair_correlation_bool):
    """
    Molecular dynamics simulation using the Euler or Verlet's algorithms
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep.

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
        
    T = 1.35
    #Create a nx3 matrix to store the velocity of each particle at each step in time.
    next_step_velocity = np.copy(init_velocity(T, n))

    #Create a nx3 matrix to store the position of each particle at each step in time.
    next_step_position = np.copy(init_pos)
    
    final_matrix_pos = np.copy(init_pos)
    final_matrix_vel = np.copy(init_velocity(T, n))
    
    rel_dist = atomic_distances(init_pos, 0)
    final_vector_kin = np.array([kin_en(init_velocity(T, n))])
    final_vector_pot = np.array([pot_en(rel_dist)])
    final_vector_energy = np.array([kin_en(init_velocity(T, n)) + pot_en(rel_dist)])
    
    final_vector_corr = pair_correlation(rel_dist, dr)
    
    if pressure_bool == True:
        final_vector_tba = np.array([pressure(atomic_distances(final_matrix_pos[:,-3:], 0))])
        final_vector_press = np.copy([final_vector_tba[0]])
        rho = n/(L*L*L)
    
    #print("Init energy:" , final_vector_energy)

    final_rel_dist = atomic_distances(init_pos, 0)
    j = 0

    for i in range(number_of_steps):
        rel_pos = atomic_distances(final_matrix_pos[:,-3:], 1)
        rel_dist = atomic_distances(final_matrix_pos[:,-3:], 0)
        next_step_position, next_step_velocity = algorithm(final_matrix_pos[:,-3:], final_matrix_vel[:,-3:], rel_pos, rel_dist)
        next_step_position = next_step_position % L
        final_matrix_pos =  np.concatenate((final_matrix_pos, next_step_position), axis=1, out=None)
        final_matrix_vel =  np.concatenate((final_matrix_vel, next_step_velocity), axis=1, out=None)
        final_vector_kin = np.concatenate((final_vector_kin, np.array([kin_en(next_step_velocity)])), axis = 0, out = None)
        final_vector_pot = np.concatenate((final_vector_pot, np.array([pot_en(rel_dist)])), axis=0, out=None)
        final_vector_energy = np.concatenate((final_vector_energy, np.array([kin_en(next_step_velocity) + pot_en(rel_dist)])), axis=0, out=None)
        final_rel_dist = np.concatenate((final_rel_dist, atomic_distances(next_step_position, 0)), axis = 1, out = None)
        if pressure_bool == True:
            T = 2*kin_en(final_matrix_vel[:,-3:])/(3*(n-1)*119.8)
            final_vector_tba = np.concatenate((final_vector_tba, np.array([pressure(rel_dist)])), axis = 0, out = None)
            final_vector_press = np.concatenate((final_vector_press, np.array([np.sum([rho*119.8/T , -(rho/(3*n))*final_vector_tba[i]])])))
        
        if pair_correlation_bool:   # will add condition that it only computes this when i>n_0 or something
            final_vector_corr += pair_correlation(rel_dist, dr)
        
        #Rescaling:
        if rescaling_bool ==True:
            window = 200
            if  i>0 and i<int(0.7*number_of_steps) and i%window == 0 and abs(final_vector_energy[-1] - np.sum(final_vector_energy[-10:])/10) > 0.01*(np.sum(final_vector_energy[-10:])/10):
                j = j+1
                l = np.sqrt((3*(n-1)*T)/(2*final_vector_kin[-1]*119.8))
                final_matrix_vel = final_matrix_vel * l
                #T = 2*kin_en(final_matrix_vel[:,-3:])/(3*(n-1)*119.8)
                print("T = ", T)
    print("Number of rescalings: ", j)
    
    # Compute pressure:
    if pressure_bool == True:
        n_0 = int(0.7*number_of_steps)
        #P = (1/(n-n_0))*np.sum(final_vector_tba[-n_0:])
        Press = np.sum([1 , -(T/119.8*(3*n))*(1/(number_of_steps-n_0))*np.sum(final_vector_tba[-n_0:])])
        print("Pressure=", Press)
        
    if pair_correlation_bool:
        final_vector_corr = final_vector_corr / number_of_steps     # number_of_steps should become n_0

    #print("Positions:\n" , final_matrix_pos)
    #print("Velocities:\n" , init_vel)
    #print("Energy:\n" , final_vector_pot)
    #print(final_vector_energy)
    #print(final_vector_pot)
    #print(final_rel_dist)

    
    #now we compute the error with the steps below
    if error_bool == True:
        N = int(0.3*number_of_steps)
        P = np.copy((final_vector_press[-N:]))
        S_a = np.zeros((int(N/6), 1))
        for b in range(1, int(N/6)):
            Nb = int(N/b)
            p=np.zeros(Nb)
            P_shape=np.reshape(P, (Nb,b))
            p = np.sum(P_shape, axis = 1)/b
            S_a[b]=np.sqrt((1/(Nb-1))*(((1/Nb)*np.sum(np.square(p)))-np.square((1/Nb)*np.sum(p))))

        
        
    return final_vector_energy, final_vector_kin, final_vector_pot, final_vector_press, final_vector_corr, Press, S_a

def plot_energy(total, kinetic, potential):
    plt.plot(times, total, label = "Total")
    plt.plot(times, kinetic, label = "Kinetic")
    plt.plot(times, potential, label = "Potential")
    plt.legend()
    plt.title("Energy over time")
    plt.xlabel("$t/(m\sigma^2/\epsilon)^{1/2}$")
    plt.ylabel("$E/\epsilon$")
    plt.show()
    
    
    #now we compute the error with the steps below
def plot_error(S_a):

    x = np.arange(0, len(S_a))
    plt.plot(x, S_a)
    plt.title("Error with data blocking")
    plt.xlabel("b")
    plt.ylabel("$\sigma_A$")
    plt.show()

def plot_pair_correlation(vector_corr):
    bins_list = np.arange(2*dr, L+2*dr, dr)
    plt.title("pair correlation function")
    plt.xlabel("r")
    plt.ylabel("g(r)")
    vector_corr = vector_corr[1:]
    x = bins_list[:(len(vector_corr))] + dr/2
    vector_corr = (2 * L**3 * vector_corr)/(4*np.pi* dr* n*(n-1)*x**2)
    plt.plot(x, vector_corr)
    plt.show()
