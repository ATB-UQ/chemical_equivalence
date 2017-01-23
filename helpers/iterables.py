from typing import List, TypeVar

from functools import reduce

T = TypeVar('T')

def concat(list_of_lists: List[List[T]]) -> List[T]:
    return reduce(
        lambda acc, e: acc + e,
        list_of_lists,
        [],
    )
