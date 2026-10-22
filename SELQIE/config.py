from dataclasses import dataclass
import numpy as np  # type: ignore[reportMissingImports]
#SqQCchem.config
#Defines objects that store config options
#Things have defaults


def hello():
    print("Hello World!")

#object for holding computation settings
@dataclass
class settings:
    method: str
    n_up: int
    n_down: int
    sel_method: str = None
    test: "str" = None
    conv: float = 1e-8
    cmin: float = 1e-3
    hci_epsilon: float = None

#object for holding result data
@dataclass
class state:
    energy: float
    n_basis: int
    basis: np.array
    civec: np.array
