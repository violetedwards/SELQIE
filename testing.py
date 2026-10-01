import SeQChem as sqc
import pyscf

sqc.config.hello()

mol = pyscf.gto.Mole()
mol.verbose = 0
mol.atom = [
    ['O', ( 0., 0.    , 0.   )],
    ['H', ( 0., -0.757, 0.587)],
    ['H', ( 0., 0.757 , 0.587)],]
mol.basis = 'sto3g'
mol.build()

settings = sqc.config.settings(n_up=2,n_down=2)

sqc.methods.selected_ci.solve(mol,settings)