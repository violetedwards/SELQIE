import SELQIE as slq
import numpy as np
import pyscf

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
settings = slq.config.settings(n_up=n_alpha,n_down=n_beta)


energy = "energy:"
print(f"my cisd {energy:<15} {slq.methods.selected_ci.solve(mol,settings)[0]:.4f} Ha")

mf = mol.RHF().run()
myci = mf.CISD().run()



print(f"PySCF cisd {energy:<12} {myci.e_tot:.4f} Ha")

