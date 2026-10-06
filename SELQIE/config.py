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
    sel_method: str
    n_up: int
    n_down: int
    test: "str" = None
    conv: float = 1e-8
    cmin: float = 1e-4

@dataclass
class selection_params:
    ACI_sigma: float = None

#object for holding result data
@dataclass
class state:
    basis: np.array
    civec: np.array
    energy: float