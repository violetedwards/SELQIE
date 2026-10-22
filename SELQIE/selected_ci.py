import os
import sys

import SELQIE as slq
import pyscf # type: ignore[reportMissingImports]
import numpy as np # type: ignore[reportMissingImports]
import scipy as sp # type: ignore[reportMissingImports]

def solve(
    mol: pyscf.gto.Mole,
    settings: slq.config.settings
) -> slq.config.state:
    """
    Solve using Selected CI.
    """
    if settings.method == "sci":
        print(f"Solving {settings.method} {settings.sel_method}")
    else:
        print(f"Solving {settings.method}")

    #Run initial HF calculation
    mf = mol.RHF().run()

    #Get system data from HF calc:
    mo_coeff = np.asarray(mf.mo_coeff)

    #get 1 and 2 electron spatial integrals
    eint1 = np.linalg.multi_dot([mo_coeff.T,np.asarray(mf.get_hcore()),mo_coeff])
    eint2 = np.asarray(pyscf.ao2mo.restore(1, pyscf.ao2mo.kernel(mf.mol, mo_coeff), mo_coeff.shape[1]))

    #get other parameters
    n_orb = eint1.shape[0]
    nelec = int(sum(mf.mol.nelec))
    spin = mf.mol.spin
    nuc_replusion = mf.mol.energy_nuc()

    #transform spatial orbital integrals into spin orbital integrals
    h_int, g_int = slq.utilities.space_to_spin_eints(eint1,eint2)

    #move to fortran for solving 
    energy, n_basis = slq.solve(settings.method, settings.sel_method, n_orb, settings.n_up, settings.n_down, h_int, g_int, nuc_replusion, settings.conv, settings.hci_epsilon or 0.0, settings.cmin)
    basis = slq._selqie.selqie_api.stored_basis
    coeffs = slq._selqie.selqie_api.stored_coeffs
    state_out = slq.config.state(energy=energy, n_basis=n_basis, basis=basis, civec=coeffs)

    #return
    return state_out

   