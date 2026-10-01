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
    compare = compare_dets(column_det,row_det)
    added = compare[0]
    removed = compare[1]
    order = len(added)

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
        return h[alpha, a] + np.sum((g[a,i,alpha,i]-g[a,i,i,alpha]) for i in row_det.occupations)
    elif order == 2:
        #dets are double excitation difference
        alpha = added[0]
        beta = added[1]
        a = removed[0]
        b = removed[1]
        return g[a,b,alpha,beta] - g[a,b,beta,alpha]

def compare_dets(det1: slq.determinants.determinant, det2: slq.determinants.determinant) -> np.array:
    """
    Compare determinants to get added/removed orbitals

    Parameters
    ----------
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