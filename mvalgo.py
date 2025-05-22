# from monord import *
from typing import Callable
from mvpoly import MvPolynomial, long_div_ls, S_pair

MO_type = Callable[[tuple[int, ...], tuple[int, ...]], int]

def check_consistency(mvp_ls: list[MvPolynomial]):
    if not mvp_ls:
        raise ValueError("The generating set is empty")
    for i in range(1, len(mvp_ls)):
        mvp_ls[0].check_type_match(mvp_ls[i])

def buchberger(basis: list[MvPolynomial], monord: MO_type):
    check_consistency(basis)
    gbasis = [f for f in basis]
    zero_poly = MvPolynomial(basis[0].varnum, basis[0].field_cls, {})
    while True:
        gbasis_cache = [f for f in gbasis]
        for i in range(len(gbasis_cache)):
            for j in range(i + 1, len(gbasis_cache)):
                si, sj = S_pair(gbasis[i], gbasis[j], monord)
                r = long_div_ls(si - sj, gbasis_cache, monord)
                if r != zero_poly:
                    gbasis.append(r)
        if len(gbasis) == len(gbasis_cache):
            break
    return gbasis

def faugere_f4(basis: list[MvPolynomial], monord: MO_type):
    check_consistency(basis)
    return NotImplementedError("TODO implement F4 algorithm")
