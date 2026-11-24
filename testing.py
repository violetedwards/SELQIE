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
C   0.000000   0.000000   1.130000
"""
mol.basis = 'sto-3g'
mol.symmetry = "c2v"
mol.build()

n_alpha, n_beta = mol.nelec
hci_settings = slq.config.settings(method="sci",sel_method="hci",n_up=n_alpha,n_down=n_beta,conv=1e-6,cmin=1e-6, hci_epsilon=0.001)
aci_settings = slq.config.settings(method="sci",sel_method="aci",n_up=n_alpha,n_down=n_beta,conv=1e-6,cmin=1e-6, aci_sigma=0.0002)
asci_settings = slq.config.settings(method="sci",sel_method="asci",n_up=n_alpha,n_down=n_beta,conv=1e-6,cmin=1e-6, asci_cdets=100, asci_tdets=575)
cipsi_settings = slq.config.settings(method="sci",sel_method="cipsi",n_up=n_alpha,n_down=n_beta,conv=1e-6,cmin=1e-6, cipsi_coeff_threshold=0.0007)


def timed(fn):
    """Run fn once, returning its result and the wall time."""
    t = time.perf_counter()
    result = fn()
    return result, time.perf_counter() - t

# every timing below is end-to-end from `mol` (HF included), so the methods are comparable
def pyscf_fci():
    mf = mol.RHF().run()
    return mf, pyscf.fci.FCI(mf).kernel()[0]

hci_state, t_slq_hci = timed(
    lambda: slq.selected_ci.solve(mol, hci_settings)
)

aci_state, t_slq_aci = timed(
    lambda: slq.selected_ci.solve(mol, aci_settings)
)

asci_state, t_slq_asci = timed(
    lambda: slq.selected_ci.solve(mol, asci_settings)
)

cipsi_state, t_slq_cipsi = timed(
    lambda: slq.selected_ci.solve(mol, cipsi_settings)
)


(mf, fci_e), t_pyscf_fci = timed(pyscf_fci)

n_orb = mf.mo_coeff.shape[1]
n_alpha, n_beta = mf.mol.nelec
raw_fci_dets = math.comb(n_orb, n_alpha) * math.comb(n_orb, n_beta)


energystr = "energy:"
timestr = "time:"
detstr = "dets:"
colon = ":"
print(f"PYSCF FCI {detstr:<11} {raw_fci_dets}")
print(f"SELQIE HCI {detstr:<10} {hci_state.n_basis}")
print(f"SELQIE ACI {detstr:<10} {aci_state.n_basis}")
print(f"SELQIE ASCI {detstr:<9} {asci_state.n_basis}")
print(f"SELQIE CIPSI {detstr:<8} {cipsi_state.n_basis}\n")

print(f"SELQIE conv {energystr:<9} {hci_settings.conv:.1g} Ha")
print(f"SELQIE HCI epsilon{colon:<3} {hci_settings.hci_epsilon} Ha")
print(f"SELQIE ACI sigma{colon:<5} {aci_settings.aci_sigma} Ha")
print(f"SELQIE ASCI cdets{colon:<4} {asci_settings.asci_cdets} dets")
print(f"SELQIE ASCI tdets{colon:<4} {asci_settings.asci_tdets} dets")
print(f"SELQIE CIPSI cutoff{colon:<2} {cipsi_settings.cipsi_coeff_threshold}\n")

print(f"PySCF FCI {energystr:<10} {fci_e:.12f} Ha")
print(f"SELQIE HCI {energystr:<9} {hci_state.energy:.12f} Ha (delta = {fci_e-hci_state.energy:.2g} Ha)")
print(f"SELQIE ACI {energystr:<9} {aci_state.energy:.12f} Ha (delta = {fci_e-aci_state.energy:.2g} Ha)")
print(f"SELQIE ASCI {energystr:<8} {asci_state.energy:.12f} Ha (delta = {fci_e-asci_state.energy:.2g} Ha)")
print(f"SELQIE CIPSI {energystr:<7} {cipsi_state.energy:.12f} Ha (delta = {fci_e-cipsi_state.energy:.2g} Ha)\n")

print(f"PySCF FCI {timestr:<11} {t_pyscf_fci:.3f} secs")
print(f"SELQIE HCI {timestr:<10} {t_slq_hci:.3f} secs")
print(f"SELQIE ACI {timestr:<10} {t_slq_aci:.3f} secs")
print(f"SELQIE ASCI {timestr:<9} {t_slq_asci:.3f} secs")
print(f"SELQIE CIPSI {timestr:<8} {t_slq_cipsi:.3f} secs\n")

