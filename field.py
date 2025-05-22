from fractions import Fraction
from galois import GF
from typing import Protocol, runtime_checkable#, TypeVar

# K = TypeVar("K", bound="Field")

@runtime_checkable
class Field(Protocol):
    def __add__(self, other: "Field") -> "Field": ...
    def __sub__(self, other: "Field") -> "Field": ...
    def __mul__(self, other: "Field") -> "Field": ...
    def __truediv__(self, other: "Field") -> "Field": ...
    def __neg__(self) -> "Field": ...
    def __eq__(self, other: object) -> bool: ...
    def __str__(self) -> str: ...

    @classmethod
    def zero(cls) -> "Field": ...
    
    @classmethod
    def one(cls) -> "Field": ...

if __name__ == "__main__":
    # patch for Fraction
    QQ = Fraction
    setattr(Fraction, 'zero', classmethod(lambda _: QQ(0)))
    setattr(Fraction, 'one', classmethod(lambda _: QQ(1)))
    GF243 = GF(3**5)
    setattr(GF243, 'zero', classmethod(lambda _: GF243(0)))
    setattr(GF243, 'one', classmethod(lambda _: GF243(1)))
    print(isinstance(QQ, Field)) # True
    print(isinstance(GF243, Field)) # True
    # Now we can do `Field.register(QQ)` and `Field.register(GF243)`