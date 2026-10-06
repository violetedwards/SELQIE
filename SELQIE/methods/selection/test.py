import numpy as np # type: ignore[reportMissingImports]
import SELQIE as slq

def select(selected_space, interacting_space, h, g, selected_c):
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
