import SELQIE as slq
import pyscf
import numpy as np

def solve(
    mol: pyscf.gto.Mole,
    settings: slq.config.settings
) -> slq.config.state:
    """
    Solve using Selected CI.
    """
    print("Solving")

    #Run initial HF calculation
    mf = pyscf.scf.RHF(mol)
    mf.kernel()

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

    #generate reference determinant
    ref_det = slq.determinants.reference_determinant(settings.n_up, settings.n_down, n_orb)

    #set initial selected space
    #get selected space hamiltonian
    #diagonalise selected hamiltonian
    #start main loop
        #generate interacting space for selected space
        #apply selection algorithm to add to seleced space
        #prune
        #get hamiltonian for new slected space
        #diagonalise new selected hamiltonian
        #check convergence
    #end main loop

    #apply PT2

    #output


