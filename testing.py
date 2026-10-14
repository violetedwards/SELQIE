import math
import time
import SELQIE as slq
import numpy as np # type: ignore[reportMissingImports]
import pyscf # type: ignore[reportMissingImports]

slq.config.hello()


mol = pyscf.gto.Mole()
mol.verbose = 0
mol.atom = """
O   0.000000   0.000000   0.000000
H   0.000000   0.758602   0.504284
H   0.000000  -0.758602   0.504284
"""
mol.basis = '3-21g'
mol.symmetry = "c2v"
mol.build()

n_alpha, n_beta = mol.nelec
settings = slq.config.settings(sel_method="hci",n_up=n_alpha,n_down=n_beta,conv=1e-6,cmin=1e-6, hci_epsilon=0.1)






t1 = time.time()
sci_e, sci_len = slq.methods.selected_ci.solve(mol,settings)
t2 = time.time()




t3 = time.time()
mf = mol.RHF().run()
cisd_obj = mf.CISD()
cisd_obj.run()
t4 = time.time()





n_orb = mf.mo_coeff.shape[1]
n_alpha, n_beta = mf.mol.nelec
raw_fci_dets = int(math.comb(n_orb, n_alpha) * math.comb(n_orb, n_beta))

t5 = time.time()
fci_obj = pyscf.fci.FCI(mf)
fci_e, _ = fci_obj.kernel()
t6 = time.time()


energystr = "energy:"
timestr = "time:"
detstr = "dets:"
print(f"\nSELQIE HCI {detstr:<10} {sci_len}")
print(f"PYSCF CISD {detstr:<10} {len(cisd_obj.ci)}")
print(f"PYSCF FCI {detstr:<11} {raw_fci_dets}\n")


print(f"SELQIE conv {energystr:<10} {settings.conv:.9f} Ha")
print(f"SELQIE HCI {energystr:<9} {sci_e:.9f} Ha")
print(f"PySCF CISD {energystr:<9} {cisd_obj.e_tot} Ha")
print(f"PySCF FCI {energystr:<10} {fci_e} Ha\n")

print(f"SELQIE HCI {timestr:<10} {t2-t1:.3f} secs")
print(f"PySCF CISD {timestr:<10} {t4-t3:.3f} secs")
print(f"PySCF FCI {timestr:<11} {t6-t5:.3f} secs\n")



