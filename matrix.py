from field import Field

def rref(field_cls: type[Field], mat: list[list[Field]]):
    rown = len(mat)
    coln = len(mat[0])
    for r in range(1, rown):
        if len(mat[r]) != coln:
            raise ValueError(f"Inconsistent row size; \
                             row 0: {coln}, row {r}: {len(mat[r])}")
    emat = [row[:] for row in mat] # deep copy of mat
    # downward
    pc, pr = 0, 0
    while pc < coln and pr < rown:
        zr = pr
        while zr < rown and emat[zr][pc] == field_cls.zero():
            zr += 1
        if zr == rown: # the column is all zero
            pc += 1
            continue
        elif zr != pr: # swap the two rows pr & zr
            tmp = emat[pr]
            emat[pr] = emat[zr]
            emat[zr] = tmp
        print(f"pivot row-col ind: {pr}, {pc}")
        # now the upper-leftmost entry of the submatrix is nonzero
        emat[pr] = [e / emat[pr][pc] for e in emat[pr]]
        for r in range(pr + 1, rown):
            emat[r] = [e - emat[r][pc] * ep for e, ep in zip(emat[r], emat[pr])]
        print("after downward elim:")
        print(emat)
        pr += 1
        pc += 1
    # upward
    for r in reversed(range(pr)):
        pc = min([c for c in range(coln) if emat[r][c] != field_cls.zero()])
        for ur in range(r):
            emat[ur] = [e - emat[ur][pc] * ep for e, ep in zip(emat[ur], emat[r])]
        print(f"after upward elim @ row-col ind {r}, {pc}:")
        print(emat)
    return emat

if __name__ == "__main__":
    from fractions import Fraction
    from typing import cast
    QQ = Fraction
    setattr(QQ, 'zero', classmethod(lambda _: QQ(0)))
    setattr(QQ, 'one', classmethod(lambda _: QQ(1)))
    QQT = cast(type[Field], QQ)
    def QQF(*args):
        return cast(Field, QQ(*args))
    
    print("[[[ mat1 rref ]]]")
    mat1 = [
        [QQF(2), QQF(1), QQF(-1), QQF(8)],
        [QQF(-3), QQF(-1), QQF(2), QQF(-11)],
        [QQF(-2), QQF(1), QQF(2), QQF(-3)],
    ]
    print(rref(QQT, mat1))
    print(mat1)

    print()
    print("[[[ mat2 rref ]]]")
    mat2 = [
        [QQF(1), QQF(0), QQF(4), QQF(2)],
        [QQF(1), QQF(2), QQF(6), QQF(2)],
        [QQF(2), QQF(0), QQF(8), QQF(8)],
        [QQF(2), QQF(1), QQF(9), QQF(4)],
    ]
    print(rref(QQT, mat2))