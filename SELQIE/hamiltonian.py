import numpy as np
import SELQIE as slq


#build hamiltonian from given ci determinant basis set
def build_hamiltonian(det_basis_set: np.array, h: np.array, g:np.array) -> np.array:
    """
    Build the hamiltonian matrix for a specific ci determinant basis set
    
    Parameters
    ----------
    det_basis_set:  np.array([determinant,...])
                    list of determinants that form the basis of the hamiltonian
    h:              np.array
                    one-electron integrals
    g:              np.array
                    two-electron integrals
                
    Returns
    -------
    hamiltonian: np.array    
            The hamiltonian in the basis set provided                   
    """
    #initialise array
    hamiltonian = np.zeros((len(det_basis_set),len(det_basis_set)),dtype=float)
    for i, column_det in enumerate(det_basis_set):  #columns
        for j in range(0, i + 1):                   #rows
            element = get_matrix_element(column_det, det_basis_set[j], h, g)
            hamiltonian[i,j] = element
            hamiltonian[j,i] = element  #hermitian so symmetric
    return hamiltonian


#calculate hamiltonain matrix element
def get_matrix_element(column_det: slq.determinants.determinant, row_det: slq.determinants.determinant, h: np.array, g: np.array) -> float:
    """
    Calculate Hamiltonian element

    Parameters
    ----------
    column_det:     slq.determinants.determinant
                    basis determimant for the column of the element being calculated
    row_det:        slq.determinants.determinant
                    basis determinant for the row of the element being calculated
    h:              np.array
                    one electron integrals
    g:              np.array
                    two electron integrals

    Returns
    -------
    element:        float
                    the calculated Hamiltonian matrix element
    """
    compare = compare_dets(row_det, column_det)
    added = compare[0]
    removed = compare[1]
    order = len(added)
    phase = excitation_phase(row_det, added, removed)

    #Slater-Condon rules
    if order > 2:
        return 0.0
    elif order == 0:
        #dets are same 
        return np.sum(h[i,i] for i in row_det.occupations) + 0.5 * np.sum((g[i,j,i,j]-g[i,j,j,i]) for i in row_det.occupations for j in row_det.occupations if i != j)
    elif order == 1:
        #dets are single exictation difference
        alpha = added[0]
        a = removed[0]
        return phase * (h[alpha, a] + np.sum((g[alpha,i,a,i]-g[alpha,i,i,a]) for i in row_det.occupations if i != a))
    elif order == 2:
        #dets are double excitation difference
        alpha = added[0]
        beta = added[1]
        a = removed[0]
        b = removed[1]
        return phase * (g[alpha,beta,a,b] - g[alpha,beta,b,a])


def excitation_phase(det: slq.determinants.determinant, added: np.array, removed: np.array) -> int:
    bitstring = det.bitstring
    phase = 1
    for orbital in removed:
        if (bitstring & ((1 << orbital) - 1)).bit_count() % 2:
            phase = -phase
        bitstring ^= 1 << orbital
    for orbital in reversed(added):
        if (bitstring & ((1 << orbital) - 1)).bit_count() % 2:
            phase = -phase
        bitstring |= 1 << orbital
    return phase

def compare_dets(det1: slq.determinants.determinant, det2: slq.determinants.determinant) -> np.array:
    """
    Compare determinants to get added/removed orbitals


    det1:           slq.determinants.determinant
                    reference determinant
    det2:           slq.determinants.determinant
                    excited determinant
    
    Returns
    -------
    orbital_diff:   np.array([added,removed])
                    added: list of orbitals addded to det2
                    removed: list of orbitals removed from det2
    """
    removed = np.array([i for i in det1.occupations if (det2.bitstring & (1<<i)) == False])
    added = np.array([i for i in det2.occupations if (det1.bitstring & (1<<i)) == False])
    return np.array([added,removed])

def space_to_spin_eints(eint1: np.array, eint2: np.array) -> tuple[np.array,np.array]:
    """
    Transform spatial 1 and 2 electron integrals so that they are in terms of spin orbitals.

    Parameters
    ----------
    eint1:      np.array
                Spacial one-electron integrals
    eint2:      np.array
                Spacial two-electron integrals

    Returns
    -------
    h_int:      np.array
                Spin-orbit one-electron integrals
    g_int:      np.array
                Spin-orbit two-electron integrals
    """

    #initialise variables
    n_orb = len(eint1)
    h_int = np.zeros((2*n_orb, 2*n_orb), dtype=float)
    g_int = np.zeros((2*n_orb, 2*n_orb, 2*n_orb, 2*n_orb), dtype=float)

    #assign h values
    for spin in range(2):
        index = slice(spin*n_orb, (spin+1)*n_orb)
        h_int[index,index] = eint1

    #assign g values
    for p in range(2*n_orb):
        for q in range(2*n_orb):
            for r in range(2*n_orb):
                for s in range(2*n_orb):
                    if (p // n_orb) == (r // n_orb) and (q //n_orb) == (s // n_orb):
                        #assign in different order. PySCF uses chemist's p,r,q,s. We will use physicist's p,q,r,s
                        g_int[p,q,r,s] = eint2[p%n_orb, r%n_orb, q%n_orb, s%n_orb]

    return h_int, g_int



