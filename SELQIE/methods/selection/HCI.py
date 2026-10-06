import numpy as np # type: ignore[reportMissingImports]
import SELQIE as slq

def select(selected_space,h,g,epsilon):
    """
    Perform HCI selection
    """

    interacting_space = slq.determinants.unique_dets(np.append(slq.determinants.get_connected_dets(selected_space,1,multiple=True),slq.determinants.get_connected_dets(selected_space,2,multiple=True)))

    new_dets = np.array([])
    for interacting_det in interacting_space:
        accept = False
        for selected_det_index in range(0,len(selected_space)):
            if (np.abs(slq.hamiltonian.get_matrix_element(interacting_det,selected_space[selected_det_index],h,g))) > epsilon:
                accept = True
        if accept == True:
            new_dets = np.append(new_dets, interacting_det)

    return slq.determinants.unique_dets(np.append(selected_space, new_dets), existing=None)