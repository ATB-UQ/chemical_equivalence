from typing import Any, Dict, Callable, Optional, Tuple, List
from logging import Logger

Atom = Dict[str, Any]

Ring = Dict[str, Any]

Coordinate = List[float]

class FlavourCounter(object):
    def __init__(self, init: int = 0) -> None:
        self.i = init
    def getNext(self) -> int:
        self.i += 1
        return self.i

MolData = Any

Exception_Searching_Function = Callable[[MolData, FlavourCounter, Logger], bool]
