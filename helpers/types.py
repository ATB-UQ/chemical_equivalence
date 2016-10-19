from typing import Any, Dict, Callable, Optional
from logging import Logger

from chemical_equivalence.molData import MolData

Atom = Dict[str, Any]

class FlavourCounter(object):
    def __init__(self, init: int = 0) -> None:
        self.i = init
    def getNext(self) -> int:
        self.i += 1
        return self.i

Exception_Searching_Function = Callable[[MolData, FlavourCounter, Logger], bool]
