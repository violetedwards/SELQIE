import SeQChem as sqc
import numpy as np
from dataclasses import dataclass



#determinant definition
@dataclass
class determinant:
    """
    stored as occupancy bitstring, 1 being occupied, 0 being unoccupied. alpha spin orbitals first, beta spin orbitals second.
    """
    bitstring: int
    n_spin_orbitals: int

    @classmethod
    def create(self, occupations: np.array, n_spin_orbitals: int) -> "determinant":
        """
        create determinant from occupation list
        """
        bitstring = 0
        for orbital in occupations:
            bitstring = bitstring | (1<<orbital) 
        return self(bitstring, n_spin_orbitals)

    def occupations(self) -> np.array:
        return np.array([i for i in range(0,self.n_spin_orbitals) if self.bitstring & (1<<i)])




#create reference determinant
def reference_determinant(n_orb: int, n_alpha: int, n_beta: int) -> determinant:
    occupations = np.concatenate((np.arange(0,n_alpha), np.arange(n_orb, n_orb+n_beta)))
    return determinant.create(occupations, 2*n_orb)