# SELQIE
SELected Quantum Interacting Electrons

# WORK IN PROGRESS - CURRENTLY NOT WORKING

A modular python package for ab initio selected quantum chemistry methods.
ReadTheDocs: https://selqie.readthedocs.io/

Install via
``pip install git+https://github.com/violetedwards/SELQIE`` 

See tutorial.ipynb for basic functionality explaination. 
See testing.py for a basic comparison to PySCF.

Requires:
Numpy
PySCF
Scipy
LAPACK
BLAS

Current code features:
- it works
- HCI

Next planned functionality (in order):
- Other Selected CI methods implemented. (ACI, ASCI, etc)
- Ability to run FORTRAN code in parallel.

Future ideas (not in order):
- Orbital optimisation
- More observables, including electromagnetic rubbish that I dont understand.
- Different initial calculations (CAS, etc)
- Selected Coupled Cluster
