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
#mol.basis = 'sto-3g'
mol.basis = "6-31g*"
mol.symmetry = "c2v"
mol.build()

n_alpha, n_beta = mol.nelec
hci_settings = slq.config.settings(method="sci",sel_method="hci",n_up=n_alpha,n_down=n_beta,conv=1e-6,cmin=1e-6, hci_epsilon=0.001)


def timed(fn):
    """Run fn once, returning its result and the wall time."""
    t = time.perf_counter()
    result = fn()
    return result, time.perf_counter() - t

# every timing below is end-to-end from `mol` (HF included), so the methods are comparable
def pyscf_fci():
    mf = mol.RHF().run()
    return mf, pyscf.fci.FCI(mf).kernel()[0]


slq_hci_out, t_slq_hci = timed(lambda: slq.methods.selected_ci.solve(mol, hci_settings)[0])

(mf, fci_e), t_pyscf_fci = timed(pyscf_fci)

n_orb = mf.mo_coeff.shape[1]
n_alpha, n_beta = mf.mol.nelec
raw_fci_dets = math.comb(n_orb, n_alpha) * math.comb(n_orb, n_beta)


energystr = "energy:"
timestr = "time:"
detstr = "dets:"
colon = ":"
print(f"\nPYSCF FCI {detstr:<11} {raw_fci_dets}")
print(f"SELQIE HCI {detstr:<10} {slq_hci_out[1]}")

print(f"SELQIE conv {energystr:<10} {hci_settings.conv:.1g} Ha")
print(f"SELQIE HCI epsilon {colon:<3} {hci_settings.hci_epsilon} Ha\n")

print(f"PySCF FCI {energystr:<10} {fci_e} Ha")
print(f"SELQIE HCI {energystr:<9} {slq_hci_out[0]:.9f} Ha\n")

print(f"PySCF FCI {timestr:<11} {t_pyscf_fci:.3f} secs")
print(f"SELQIE HCI {timestr:<10} {t_slq_hci:.3f} secs\n")



