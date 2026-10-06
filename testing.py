import math
import SELQIE as slq
import numpy as np # type: ignore[reportMissingImports]
import pyscf # type: ignore[reportMissingImports]

slq.config.hello()

mol = pyscf.gto.Mole()
mol.verbose = 0
mol.atom = [
    ['O', ( 0., 0.    , 0.   )],
    ['H', ( 0., -0.757, 0.587)],
    ['H', ( 0., 0.757 , 0.587)],]
mol.basis = 'sto3g'
mol.build()

n_alpha, n_beta = mol.nelec
settings = slq.config.settings(sel_method="hci",n_up=n_alpha,n_down=n_beta,conv=1e-9,cmin=1e-6, hci_epsilon=0.00001)

energy = "energy:"

sci_e, sci_len = slq.methods.selected_ci.solve(mol,settings)
print(f"SELQIE hci dets: {sci_len}")
print(f"SELQIE converged to {energy:<8} {settings.conv:.9f} Ha")
print(f"SELQIE hci {energy:<15} {sci_e} Ha\n")

mf = mol.RHF().run()

cisd_obj = mf.CISD()
cisd_obj.run()
print(f"PYSCF cisd dets: {len(cisd_obj.ci)}")
print(f"PySCF cisd {energy:<15} {cisd_obj.e_tot} Ha\n")

n_orb = mf.mo_coeff.shape[1]
n_alpha, n_beta = mf.mol.nelec
raw_fci_dets = int(math.comb(n_orb, n_alpha) * math.comb(n_orb, n_beta))
print(f"PYSCF fci dets: {raw_fci_dets}")

fci_obj = pyscf.fci.FCI(mf)
fci_e, _ = fci_obj.kernel()
print(f"PySCF fci {energy:<16} {fci_e} Ha")
