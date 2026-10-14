import os
import sys

import SELQIE as slq
import pyscf # type: ignore[reportMissingImports]
import numpy as np # type: ignore[reportMissingImports]
import scipy as sp # type: ignore[reportMissingImports]

_fortran_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'Fortran', 'build'))
if os.path.isdir(_fortran_dir) and _fortran_dir not in sys.path:
    sys.path.insert(0, _fortran_dir)

try:
    _selected_ci = __import__('selected_ci')
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "The Fortran extension is not built yet. Run 'meson setup build && meson compile -C build' in the Fortran directory first."
    ) from exc

def add_selected_dets(algorithm: str,selected_space: np.array,selected_E: float,selected_c: np.array, h: np.array, g: np.array, settings:slq.config.settings) -> np.array:
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
    if algorithm == "test":
        return slq.methods.selection.test.select(selected_space,h,g,selected_c)
    if algorithm == "hci":
        return slq.methods.selection.HCI.select(selected_space,h,g,settings.hci_epsilon)



def solve(
    mol: pyscf.gto.Mole,
    settings: slq.config.settings
) -> slq.config.state:
    """
    Solve using Selected CI.
    """
    print("Solving")

    #Run initial HF calculation
    mf = mol.RHF().run()

    #Get system data from this:
    mo_coeff = np.asarray(mf.mo_coeff)
    #1 and 2 electron integrals
    eint1 = np.linalg.multi_dot([mo_coeff.T,np.asarray(mf.get_hcore()),mo_coeff])
    eint2 = np.asarray(pyscf.ao2mo.restore(1, pyscf.ao2mo.kernel(mf.mol, mo_coeff), mo_coeff.shape[1]))
    #other parameters
    n_orb = eint1.shape[0]
    nelec = int(sum(mf.mol.nelec))
    spin = mf.mol.spin
    nuc_replusion = mf.mol.energy_nuc()

    #transform spatial orbital integrals into spin orbital integrals
    h_int, g_int = slq.hamiltonian.space_to_spin_eints(eint1,eint2)

    energy = _selected_ci.solve(settings.method, settings.sel_method, n_orb, settings.n_up, settings.n_down, h_int, g_int, nuc_replusion)
    return energy, 1

    #----------------------------------------

    #generate reference determinant
    ref_det = slq.determinants.reference_determinant(n_orb, settings.n_up, settings.n_down)

  

    #set initial selected space
    selected_space = slq.determinants.initial_selected_space(ref_det)
    print(f"Initial Selected space size: {len(selected_space)}")
    #start main loop
    converged = False
    iteration = 0
    while converged == False:
        print(f"Iteraion {iteration}")
        #get selected space hamiltonian
        selected_space_H = slq.hamiltonian.build_hamiltonian(selected_space, h_int, g_int)

        #diagonalise selected hamiltonian
        selected_E_old, selected_c_old = sp.linalg.eigh(selected_space_H)
        if type(selected_E_old) == np.ndarray:
            selected_E_old = selected_E_old[0]
            selected_c_old = selected_c_old[:,0]

        #apply selection algorithm to add to seleced space
        #print("performing selection...")
        selected_space = add_selected_dets(settings.sel_method, selected_space, selected_E_old, selected_c_old, h_int, g_int,settings)

        #get hamiltonian for new slected space
        selected_space_H = slq.hamiltonian.build_hamiltonian(selected_space, h_int, g_int)

        #diagonalise new selected hamiltonian
        selected_E_new, selected_c_new = sp.linalg.eigh(selected_space_H)
        if type(selected_E_new) == np.ndarray:
            selected_E_new = selected_E_new[0]
            selected_c_new  =selected_c_new[:,0]

        print(f"\tSelected space size: {len(selected_space)}")

        #check convergence
        print(f"\told: {selected_E_old}, new: {selected_E_new}")
        if abs(selected_E_old - selected_E_new) < settings.conv:
            converged = True
        else:
            converged = False

        #prune if not converged
        if converged == False:
            print("\tPruning...")
            keep = np.abs(selected_c_new) >= settings.cmin
            selected_space = selected_space[keep]


        print(f"\tSelected space size: {len(selected_space)}")
        iteration+=1

    #end main loop
    print("energy converged")
    print("---------------------------------------------------------------------------------\n\n")
    #apply PT2

    #output
    return np.round(selected_E_new + nuc_replusion, decimals=slq.utilities.find_exp(settings.conv)),len(selected_space)


