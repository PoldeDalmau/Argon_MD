
def lj_force(position):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    position : np.ndarray
        Relative particle positions as obtained from atomic_distances.
    

    Returns
    -------
    np.ndarray
        The net force acting on particle i due to all other particles.
    """
    rel_pos = atomic_distances(position[:,-3:], 1) # calculates the relative position of the particles, with pos = np.array(N_particels, 3)
    r = atomic_distances(position[:,-3:], 0) # n x n simmetric matrix, r[i,j] is the distance between the i-th and the j-th particles
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
