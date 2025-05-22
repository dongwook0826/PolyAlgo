from fractions import Fraction as QQ
# from galois import GF
from field import Field
from typing import Self, cast, Callable#, TypeVar, Generic, Dict, Tuple, Literal
from monord import mon_lcm, divisible_by, lex, grlex, grevlex
from functools import cmp_to_key

# VN = TypeVar('VN', bound = int)
# K = TypeVar('K', bound = Field)
MO_type = Callable[[tuple[int, ...], tuple[int, ...]], int]

class MvPolynomial: # ~~~(Generic[VN, K]): # removed generic typing parts
    # vector of exponents |-> coefficient
    def __init__(self,
                 varnum: int,
                 field_cls: type[Field],
                 poly: dict[tuple[int, ...], Field]):
        # if isinstance(poly, dict):
        poly_new = {}
        for m, c in poly.items():
            if len(m) != varnum:
                raise ValueError(f"The monomial {m} doesn\'t match varnum = {varnum}")
            if c != field_cls.zero():
                poly_new[m] = c
        self.poly: dict[tuple[int, ...], Field] = poly_new
        self.varnum = varnum
        self.field_cls = field_cls

    def type_matches(self, poly2: Self):
        return self.varnum == poly2.varnum and self.field_cls == poly2.field_cls

    def check_type_match(self, poly2: Self):
        if not self.type_matches(poly2):
            raise ValueError(f"The operand MvPolynomials\' types don\'t match")

    def __neg__(self):
        poly_neg = self.poly.copy()
        for m, c in poly_neg.items():
            poly_neg[m] = -c
        return MvPolynomial(self.varnum, self.field_cls, poly_neg)
    
    def __add__(self, poly2: Self):
        self.check_type_match(poly2)
        poly_sum = self.poly.copy()
        for m, c in poly2.poly.items():
            if m in poly_sum:
                poly_sum[m] += c
            else:
                poly_sum[m] = c
        return MvPolynomial(self.varnum, self.field_cls, poly_sum)
    
    def __sub__(self, poly2: Self):
        self.check_type_match(poly2)
        poly_diff = self.poly.copy()
        for m, c in poly2.poly.items():
            if m in poly_diff:
                poly_diff[m] -= c
            else:
                poly_diff[m] = -c
        return MvPolynomial(self.varnum, self.field_cls, poly_diff)

    def __mul__(self, poly2: Self):
        self.check_type_match(poly2)
        poly_prod = {}
        for m1, c1 in self.poly.items():
            for m2, c2 in poly2.poly.items():
                mp = tuple(x + y for x, y in zip(m1, m2))
                if mp in poly_prod:
                    poly_prod[mp] += c1 * c2
                else:
                    poly_prod[mp] = c1 * c2
        return MvPolynomial(self.varnum, self.field_cls, poly_prod)

    def __eq__(self, obj: object):
        if not isinstance(obj, MvPolynomial):
            return False
        obj = cast(Self, obj)
        for m, c in self.poly.items():
            if m not in obj.poly or obj.poly[m] != c:
                return False
        return self.type_matches(obj)
    
    def __str__(self):
        poly_str = "{"
        for i, (m, c) in enumerate(self.poly.items()):
            poly_str += (" " if i == 0 else ", ") + str(m) + ": " + str(c)
        poly_str += " }"
        return poly_str
    
    def coeff(self, mon: tuple[int, ...]):
        if len(mon) != self.varnum:
            raise ValueError(f"The monomial {mon} doesn\'t match varnum = {self.varnum}")
        return self.poly[mon] if mon in self.poly else self.field_cls.zero()
    
    def lead_mon(self, monord: MO_type):
        return max(self.poly, key = cmp_to_key(monord))
    
    def mon_list(self, monord: MO_type):
        return sorted(self.poly, key = cmp_to_key(monord))

    
# long division result (quotient & remainder) of MvPoly's
def long_div(poly1: MvPolynomial, poly2: MvPolynomial, # dividing poly1 by poly2
             monord: MO_type):
    poly1.check_type_match(poly2)
    varnum, field_cls = poly1.varnum, poly1.field_cls
    qnt = MvPolynomial(varnum, field_cls, {})
    rmd = MvPolynomial(varnum, field_cls, poly1.poly.copy())
    p2_lead_mon = poly2.lead_mon(monord)
    p2_lead_coeff = poly2.coeff(p2_lead_mon)
    while True:
        divisible_mons = [m for m in rmd.poly if divisible_by(m, p2_lead_mon)]
        if not divisible_mons:
            break
        pivot_mon = max(divisible_mons)
        qnt_mon = tuple([ep - em for ep, em in zip(pivot_mon, p2_lead_mon)])
        qnt_term = MvPolynomial(
            varnum, field_cls,{ qnt_mon: rmd.coeff(pivot_mon) / p2_lead_coeff })
        qnt += qnt_term
        rmd -= qnt_term * poly2
    return qnt, rmd

def long_div_ls(poly1: MvPolynomial, polys2: list[MvPolynomial],
                monord: MO_type):
    for poly2 in polys2:
        poly1.check_type_match(poly2)
    zero_poly = MvPolynomial(poly1.varnum, poly1.field_cls, {})
    rmd = MvPolynomial(poly1.varnum, poly1.field_cls, poly1.poly.copy())
    ps2_lm_qnt = [[poly2,
                   poly2.lead_mon(monord),
                   MvPolynomial(poly1.varnum, poly1.field_cls, {})]
                   for poly2 in polys2 if poly2 != zero_poly]
    while True:
        divided = False
        for p2, lm, qnt in ps2_lm_qnt:
            if divisible_by(rmd.lead_mon(monord), lm):
                q, r = long_div(rmd, p2, monord)
                qnt += q
                rmd = r
                divided = True
                break
        if not divided: # no more divisible poly's
            break
    return [[poly2, qnt] for poly2, _, qnt in ps2_lm_qnt], rmd
# TODO: test long_div_ls

def S_pair(poly1: MvPolynomial, poly2: MvPolynomial,
           monord: MO_type):
    poly1.check_type_match(poly2)
    lm1, lm2 = poly1.lead_mon(monord), poly2.lead_mon(monord)
    lcm = mon_lcm(lm1, lm2)
    dm1 = tuple([em - e1 for em, e1 in zip(lcm, lm1)])
    dm2 = tuple([em - e2 for em, e2 in zip(lcm, lm2)])
    s1 = poly1 * MvPolynomial(poly1.varnum, poly1.field_cls,
                              { dm1: poly1.field_cls.one() / poly1.poly[lm1] })
    s2 = poly2 * MvPolynomial(poly2.varnum, poly2.field_cls,
                              { dm2: poly2.field_cls.one() / poly2.poly[lm2] })
    return s1, s2


if __name__ == "__main__":
    
    from fractions import Fraction as QQ
    setattr(QQ, 'zero', classmethod(lambda _: QQ(0)))
    setattr(QQ, 'one', classmethod(lambda _: QQ(1)))
    # Field.register(QQ) # don't use this
    QQT = cast(type[Field], QQ)
    def QQF(*args):
        return cast(Field, QQ(*args))
    print(QQF(1))
    print(QQF(2, 3))
    print(QQF(-5))

    f = MvPolynomial(3, QQT, {
        (2,0,0): QQF(1),
        (1,1,0): QQF(1, 2),
        (0,0,0): QQF(-1),
    })
    print(f)

    f1 = MvPolynomial(3, QQT, {
        (1,1,0): QQF(1, 2),
        (0,0,0): QQF(-1),
        (2,0,0): QQF(1),
    })
    f2 = MvPolynomial(3, QQT, {
        (2,0,0): QQF(1),
        (0,0,2): QQF(-1),
    })
    f3 = MvPolynomial(3, QQT, {
        (1,1,0): QQF(1),
        (0,0,0): QQF(1),
    })
    print(f1 + f2)
    print(f1 - f2)
    print(-f3)
    print(f1 * f3)
    print(f == f1) # True
    print(f2 != f3) # True
    print(f1 == f3) # False

    g = MvPolynomial(2, QQT, {
        (1,1): QQF(1, 2),
        (0,0): QQF(-1),
        (2,0): QQF(1)
    })
    print(f == g) # False

    q1, r1 = long_div(f1, f3, grevlex)
    print(q1)
    print(r1)
    q2, r2 = long_div(f2, f3, grlex)
    print(q2)
    print(r2)