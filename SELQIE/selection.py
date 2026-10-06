import SELQIE as slq
import numpy as np
import pyscf



def add_to_selected_space(algorithm: str,selected_space: np.array,interacting_space: np.array,selected_E: float,selected_c: np.array, h: np.array, g: np.array) -> np.array:
    """
    Performs selection algorithm to decide elements of interacting space to add to selected space

    Parameters:
    algorithm           str
                        the selection algorithm to be performed
    selected_space      np.array[slq.determimants.determinant, ...]
                        the space of selected determinants before the anything is added
    interacting_space   np.array[slq.determimants.determinant, ...]
                        the space of elements that can interact with elements of the selected space
    selected_E          float
                        the energy of the solved selected space state
    selected_c          np.array[float, ...]
                        the c vector of the solved selected space state

    Returns:
    new_selected_space  np.array[slq.determimants.determinant, ...]
                        the selected space after additions
    """
    if algorithm == "dumb":
        return dumb(selected_space, interacting_space,h,g,selected_c)


def dumb(selected_space, interacting_space, h, g, selected_c):
    """
    Minimal selected-CI fix:
    keep candidates with the largest |<det | H | current_space>| weight.
    """
    if len(interacting_space) == 0:
        return selected_space

    scores = []
    for det in interacting_space:
        score = 0.0
        for i, ref_det in enumerate(selected_space):
            amp = selected_c[i]
            score += abs(slq.hamiltonian.get_matrix_element(det, ref_det, h, g) * amp)
        scores.append(score)

    order = np.argsort(scores)[::-1]
    keep = min(max(1, len(interacting_space) // 5), len(interacting_space))
    new_dets = interacting_space[order[:keep]]

    # avoid duplicates
    out = []
    seen = {d.bitstring for d in selected_space}
    for d in new_dets:
        if d.bitstring not in seen:
            out.append(d)
            seen.add(d.bitstring)

    return np.array(list(selected_space) + out, dtype=object)


