from typing import Any, Dict, Callable, Optional, Tuple, List
from logging import Logger

from chemical_equivalence.molData import MolData

Atom = Dict[str, Any]

Ring = Dict[str, Any]

Coordinate = List[float]

class FlavourCounter(object):
    def __init__(self, init: int = 0) -> None:
        self.i = init
    def getNext(self) -> int:
        self.i += 1
        return self.i

Exception_Searching_Function = Callable[[MolData, FlavourCounter, Logger], bool]
