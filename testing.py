import contextlib
import io
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
mol.basis = 'sto-3g'
mol.symmetry = "c2v"
mol.build()

n_alpha, n_beta = mol.nelec
settings = slq.config.settings(method="cisd",n_up=n_alpha,n_down=n_beta,conv=1e-9,cmin=1e-6, hci_epsilon=0.1)
fci_settings = slq.config.settings(method="fci",n_up=n_alpha,n_down=n_beta,conv=1e-9,cmin=1e-6, hci_epsilon=0.1)






def timed(fn, repeats=5):
    """Run fn repeats times, returning its result and the best wall time."""
    best = float("inf")
    for _ in range(repeats):
        t = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            result = fn()
        best = min(best, time.perf_counter() - t)
    return result, best

# every timing below is end-to-end from `mol` (HF included), so the methods are comparable
def pyscf_cisd():
    mf = mol.RHF().run()
    return mf, mf.CISD().run()

def pyscf_fci():
    mf = mol.RHF().run()
    return mf, pyscf.fci.FCI(mf).kernel()[0]

# untimed warm-up so one-time library start-up costs are not charged to the first method
slq.methods.selected_ci.solve(mol, settings)
pyscf_cisd()
pyscf_fci()

slq_out, t_slq_cisd = timed(lambda: slq.methods.selected_ci.solve(mol, settings)[0])
slq_fci_out, t_slq_fci = timed(lambda: slq.methods.selected_ci.solve(mol, fci_settings)[0])
(mf, cisd_obj), t_pyscf_cisd = timed(pyscf_cisd)
(_, fci_e), t_pyscf_fci = timed(pyscf_fci)

n_orb = mf.mo_coeff.shape[1]
n_alpha, n_beta = mf.mol.nelec
raw_fci_dets = math.comb(n_orb, n_alpha) * math.comb(n_orb, n_beta)


energystr = "energy:"
timestr = "time:"
detstr = "dets:"
print(f"\nSELQIE CISD {detstr:<9} {slq_out[1]}")
print(f"PYSCF CISD {'CSFs:':<10} {len(cisd_obj.ci)}")
print(f"SELQIE FCI {detstr:<10} {slq_fci_out[1]}")
print(f"PYSCF FCI {detstr:<11} {raw_fci_dets}\n")


print(f"SELQIE conv {energystr:<10} {settings.conv:.9f} Ha")
print(f"SELQIE CISD {energystr:<8} {slq_out[0]:.9f} Ha")
print(f"PySCF CISD {energystr:<9} {cisd_obj.e_tot} Ha")
print(f"SELQIE FCI {energystr:<9} {slq_fci_out[0]:.9f} Ha")
print(f"PySCF FCI {energystr:<10} {fci_e} Ha\n")

print(f"SELQIE CISD {timestr:<9} {t_slq_cisd:.3f} secs")
print(f"PySCF CISD {timestr:<10} {t_pyscf_cisd:.3f} secs")
print(f"SELQIE FCI {timestr:<10} {t_slq_fci:.3f} secs")
print(f"PySCF FCI {timestr:<11} {t_pyscf_fci:.3f} secs\n")



