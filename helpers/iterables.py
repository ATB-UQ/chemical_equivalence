from typing import List, Any

from functools import reduce

def concat(list_of_lists: List[List[Any]]) -> List[Any]:
    return reduce(
        lambda acc, e: acc + e,
        list_of_lists,
        [],
    )
