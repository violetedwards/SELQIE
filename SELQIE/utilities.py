from math import log10, floor
import numpy as np # type: ignore[reportMissingImports]

def find_exp(number) -> int:
    base10 = log10(abs(number))
    return abs(floor(base10))

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