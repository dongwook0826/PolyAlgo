from typing import Callable
from mvpoly import MvPolynomial, long_div_ls, S_pair
from monord import divisible_by, mon_div, mon_lcm
from functools import reduce, cmp_to_key
from matrix import rref

MO_type = Callable[[tuple[int, ...], tuple[int, ...]], int] # monomial order
choice_type = Callable[[list[MvPolynomial], set[tuple[int, int]], MO_type],
                       set[tuple[int, int]]] # pair selection strategy in F4

def check_consistency(mvp_ls: list[MvPolynomial]):
    if not mvp_ls:
        raise ValueError("The generating set is empty")
    for i in range(1, len(mvp_ls)):
        mvp_ls[0].check_type_match(mvp_ls[i])

### Buchberger algorithm ###

def buchberger(basis: list[MvPolynomial],
               monord: MO_type) -> list[MvPolynomial]:
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
        if len(gbasis) == len(gbasis_cache): # no more extension in gbasis
            break
    return gbasis

### Faugere F4 algorithm ###

def mvpoly_mat_manip(ext_s_halves: list[MvPolynomial],
                     monord: MO_type) -> list[MvPolynomial]:
    if not ext_s_halves:
        raise ValueError("The input mvpoly list to obtain RREF is empty")
    varnum = ext_s_halves[0].varnum
    field_cls = ext_s_halves[0].field_cls
    zero_poly = MvPolynomial(varnum, field_cls, {})
    # for f in ext_s_halves:
    #     if f == zero_poly:
    #         raise ValueError("A 0-polynomial is in RREF subroutine input")
    ext_s_halves = [f for f in ext_s_halves if f != zero_poly]
    zero = field_cls.zero()
    sorted_mon = sorted(reduce(lambda s, t: s | t,
                               [{m for m in p.poly} for p in ext_s_halves],
                               set()),
                        key = cmp_to_key(monord))
    lead_ind = lambda ls: min(filter(lambda i: ls[i] != zero, range(len(ls))))
    esh_mat = [[p.coeff(m) for m in sorted_mon] for p in ext_s_halves]
    lmi_esh_mat = {lead_ind(r) for r in esh_mat}
    red_esh_mat = rref(ext_s_halves[0].field_cls, esh_mat)
    red_esh = []
    for row in red_esh_mat:
        if all(map(lambda e: e == zero, row)):
            break
        if lead_ind(row) not in lmi_esh_mat:
            red_esh.append(MvPolynomial(
                varnum, field_cls,
                {m: c for m, c in zip(sorted_mon, row) if c != zero})
            )
    return red_esh

def symbolic_preprocess(s_halves: list[MvPolynomial],
                        curr_gbasis: list[MvPolynomial],
                        monord: MO_type) -> list[MvPolynomial]:
    if not s_halves:
        raise ValueError("The input mvpoly list to preprocess is empty")
    if not curr_gbasis:
        raise ValueError("The initial basis set for preprocessing is empty")
    ext_s_halves = [f for f in s_halves]
    done_mons = {f.lead_mon(monord) for f in ext_s_halves}
    while True:
        all_mons = reduce(lambda s, t: s | t,
                          [{m for m in p.poly} for p in ext_s_halves],
                          set())
        if all_mons == done_mons:
            break
        newmon = max(all_mons - done_mons, key = cmp_to_key(monord))
        done_mons |= {newmon}
        for f in curr_gbasis:
            lmf = f.lead_mon(monord)
            if divisible_by(newmon, lmf):
                ext_s_halves.append(f * MvPolynomial(
                    f.varnum, f.field_cls,
                    {mon_div(newmon, lmf): f.field_cls.one()}
                ))
    return mvpoly_mat_manip(ext_s_halves, monord)

def pair_deg(f1: MvPolynomial, f2: MvPolynomial, monord: MO_type) -> int:
    return sum(mon_lcm(f1.lead_mon(monord), f2.lead_mon(monord)))

def normal_selection(basis: list[MvPolynomial],
                     ipairs: set[tuple[int, int]],
                     monord: MO_type) -> set[tuple[int, int]]:
    ip_deg = [pair_deg(basis[i], basis[j], monord) for i, j in ipairs]
    min_deg = min(ip_deg)
    return set([(i, j) for (i, j), d in zip(ipairs, ip_deg) if d == min_deg])

def faugere_f4(basis: list[MvPolynomial],
               pair_choice: choice_type,
               monord: MO_type) -> list[MvPolynomial]:
    check_consistency(basis)
    gbasis = [f for f in basis]
    ipairs: set[tuple[int, int]] = set()
    for j in range(len(gbasis)):
        ipairs |= {(i, j) for i in range(j)}
    while ipairs:
        ip_sub = pair_choice(gbasis, ipairs, monord)
        ipairs -= ip_sub
        s_halves = reduce(
            lambda ss, s: ss + s, 
            map(lambda i: list(S_pair(gbasis[i[0]], gbasis[i[1]], monord)),
                ip_sub),
            [],
        )
        red_esh = symbolic_preprocess(s_halves, gbasis, monord)
        for k in range(len(red_esh)):
            ipairs |= {(i, len(gbasis) + k) for i in range(len(gbasis) + k)}
        gbasis += red_esh
    return gbasis

if __name__ == "__main__":
    # TODO: test the algorithms
    pass