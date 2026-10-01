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



