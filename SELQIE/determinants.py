import SELQIE as slq
import numpy as np
from dataclasses import dataclass
import itertools



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

    @property
    def occupations(self) -> np.array:
        return np.array([i for i in range(0,self.n_spin_orbitals) if self.bitstring & (1<<i)])


#create reference determinant
def reference_determinant(n_orb: int, n_alpha: int, n_beta: int) -> determinant:
    occupations = np.append(np.arange(0,n_alpha), np.arange(n_orb, n_orb+n_beta))
    return determinant.create(occupations, 2*n_orb)



def get_connected_dets(ref_det: determinant, order: int) -> np.array:
    """
    Get all determinants that can be gained by excitation (connected) from a reference by exchanging n=order electrons
    
    Parameters
    ----------
    ref_det:    slq.determinants.determinant
                determinant which will act as the reference
    order:      int
                max order of excitations from reference that will be generated.
    
    Returns
    -------
    connected_dets: np.array([slq.determinants.determinant])
                    list of connected determinants
    """
    occupied = ref_det.occupations
    virtual = np.array([i for i in range(0,ref_det.n_spin_orbitals) if (ref_det.bitstring & (1 << i))==False])
    n_orb = ref_det.n_spin_orbitals // 2
    connected_dets = np.array([])
    for removed in itertools.combinations(occupied, order):
        for added in itertools.combinations(virtual,order):
            if sum(orbital < n_orb for orbital in removed) != sum(orbital < n_orb for orbital in added):
                continue
            bitstring = ref_det.bitstring
            for orbital in removed:
                bitstring = bitstring ^ (1<<orbital)
            for orbital in added:
                bitstring = bitstring | (1<<orbital)
            connected_dets = np.append(connected_dets, determinant(bitstring, ref_det.n_spin_orbitals))

    return connected_dets