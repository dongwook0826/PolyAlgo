def check_len_match(m1: tuple[int, ...], m2: tuple[int, ...]):
    if len(m1) != len(m2):
        raise ValueError(f"Monomial sizes don't match: {len(m1)} vs {len(m2)}")
    
def mon_lcm(m1: tuple[int, ...], m2: tuple[int, ...]):
    check_len_match(m1, m2)
    return tuple([max(e1, e2) for e1, e2 in zip(m1, m2)])

def divisible_by(m1: tuple[int, ...], m2: tuple[int, ...]):
    check_len_match(m1, m2)
    return all([e1 >= e2 for e1, e2 in zip(m1, m2)])

def mon_div(m1: tuple[int, ...], m2: tuple[int, ...]):
    check_len_match(m1, m2)
    return tuple([e1 - e2 for e1, e2 in zip(m1, m2)])

### monomial orders ###
# LT = -1
# EQ = 0
# GT = 1

def lex(m1: tuple[int, ...], m2: tuple[int, ...]):
    for p in range(len(m1)):
        if m1[p] != m2[p]:
            return m1[p] - m2[p]
    return 0
# e.g. x > y > z, x^2 > xy^2 > xz > y^100 > y^2z > y^2 > z^100

def grlex(m1: tuple[int, ...], m2: tuple[int, ...]):
    gr1, gr2 = sum(m1), sum(m2)
    if gr1 != gr2:
        return gr1 - gr2
    return lex(m1, m2)
# e.g. x^2 > xy > xz > y^2 > yz > z^2 > x > y > z > 1

def grevlex(m1: tuple[int, ...], m2: tuple[int, ...]):
    gr1, gr2 = sum(m1), sum(m2)
    if gr1 != gr2:
        return gr1 - gr2
    for p in reversed(range(len(m1))):
        if m1[p] != m2[p]:
            return m2[p] - m1[p]
    return 0
# e.g. x^2 > xy > y^2 > xz > yz > z^2 > x > y > z > 1

if __name__ == "__main__":
    from functools import cmp_to_key
    for cmp in [lex, grlex, grevlex]:
        print(sorted([
            (0,0,0),
            (1,0,0), (0,1,0), (0,0,1),
            (2,0,0), (0,2,0), (0,0,2), (1,1,0), (1,0,1), (0,1,1),
            (3,0,0), (0,3,0), (0,0,3), (2,1,0), (1,2,0), (2,0,1), (1,0,2), (0,2,1), (0,1,2), (1,1,1)
        ], key = cmp_to_key(cmp), reverse = True))
    