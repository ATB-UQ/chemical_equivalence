from typing import Dict, Callable, Optional, Tuple, List
from logging import Logger

from atb_outputs.mol_data import Atom, Ring, Coordinate, MolData

class FlavourCounter(object):
    def __init__(self, init: int = 0) -> None:
        self.i = init
    def getNext(self) -> int:
        self.i += 1
        return self.i

Exception_Searching_Function = Callable[[MolData, FlavourCounter, Logger], bool]
