import SELQIE as slq
import numpy as np # type: ignore[reportMissingImports]
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



def get_connected_dets(ref_dets: determinant, order: int, multiple: bool = False) -> np.array:
    """
    Get all determinants that can be gained by excitation (connected) from a reference by exchanging n=order electrons
    
    Parameters
    ----------
    ref_det:    slq.determinants.determinant
                determinant which will act as the reference
    order:      int
                max order of excitations from reference that will be generated.
    multiple:   are there multile reference dets
    Returns
    -------
    connected_dets: np.array([slq.determinants.determinant])
                    list of connected determinants
    """
    connected_dets = np.array([])
    if multiple == True:
        for ref_det in ref_dets:
            occupied = ref_det.occupations
            virtual = np.array([i for i in range(0,ref_det.n_spin_orbitals) if (ref_det.bitstring & (1 << i))==False])
            n_orb = ref_det.n_spin_orbitals // 2
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
    if multiple == False:
        ref_det = ref_dets
        occupied = ref_det.occupations
        virtual = np.array([i for i in range(0,ref_det.n_spin_orbitals) if (ref_det.bitstring & (1 << i))==False])
        n_orb = ref_det.n_spin_orbitals // 2
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

    return unique_dets(connected_dets)

def initial_selected_space(ref_det: determinant) -> np.array:
    """
    Generate the initial selected space  ==CURRENTLY SIMPLISTIC IMPLEMENTENTATION==

    Parameters
    ----------
    ref_det     slq.determinants.determinant
                reference hartree-fock state
    
    Returns
    -------
    init_sel_space  np.array[slq.determinmants.determinant,...]
                    initial selected space
    """
    return np.array([ref_det])


def unique_dets(det_list, existing=None):
    """
    Return the determinant list with duplicates removed while preserving the
    first occurrence of each bitstring.

    ``existing`` may be passed to exclude determinants that are already in the
    selected space.
    """
    if det_list is None or len(det_list) == 0:
        return np.array([], dtype=object)

    seen = np.array([])
    if existing is not None:
        for det in existing:
            if det is not None:
                seen = np.append(seen,det.bitstring)

    unique_det_list = np.array([])
    for det in det_list:
        if det is None:
            continue
        if det.bitstring in seen:
            continue
        seen = np.append(seen,det.bitstring)
        unique_det_list = np.append(unique_det_list,det)

    return unique_det_list