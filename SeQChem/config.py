from dataclasses import dataclass
import numpy as np
#SqQCchem.config
#Defines objects that store config options
#Things have defaults


def hello():
    print("Hello World!")

#object for holding computation settings
@dataclass
class settings:
    n_up: int
    n_down: int

#object for holding result data
@dataclass
class state:
    basis: np.array
    civec: np.array
    energy: float