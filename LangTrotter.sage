def l_adic_valuation(n,l): # input: an integer n and a prime l. Output: the l-adic valuation of n
    e = 1
    while n % (l^e) == 0:
        e += 1
    return e-1


def density_LT(E,P,x): # input: an elliptic curve E, a point P of E of infinite order and a real number x>0. Output: proportion of prime numbers p<x such that P mod p generates E(F_p)
    N = E.conductor()
    n = 0
    nb_primes = 0
    for p in range(2,x):
        if is_prime(p) and N%p != 0:
            nb_primes += 1
            Ep = E.reduction(p)
            Pp = P.reduction(p)
            if Pp.order() == Ep.order():
                n+=1
    return n/nb_primes

def is_divisible(E,P,p): # returns True if there is a prime ell such that ell divides #E(F_p) and P mod p is ell-divisible, else returns False
    Ep = E.reduction(p)
    Pp = P.reduction(p)
    n = Ep.order()
    nP = Pp.order()
    if is_prime(n):
        if nP == 1:
            return True
        else:
            return False
    for l in range(2,n+1):
        if is_prime(l) and n%l ==0:
            if l_adic_valuation((Ep.abelian_group()).exponent(),l) != l_adic_valuation(nP,l):
                if Pp.division_points(l) != []:
                    return True
    return False

def is_divisible_without_2(E,P,p): # returns True if there is a prime ell > 2 such that ell divides #E(F_p) and P mod p is ell-divisible, else returns False
    Ep = E.reduction(p)
    Pp = P.reduction(p)
    n = Ep.order()
    nP = Pp.order()
    if is_prime(n):
        if nP == 1:
            return True
        else:
            return False
    for l in range(3,n+1):
        if is_prime(l) and n%l ==0:
            if l_adic_valuation((Ep.abelian_group()).exponent(),l) != l_adic_valuation(nP,l):
                if Pp.division_points(l) != []:
                    return True
    return False


def density_LT_indivisible_without_2(E,P,x): # input: an elliptic curve E, a point P of E of infinite order and a real number x>0. Output: proportion of prime numbers p<x such that the Indivisibility condition holds at all primes ell different from 2
    N = E.conductor()
    n = 0
    nb_primes = 0
    for p in range(2,x):
        if is_prime(p) and N%p != 0:
            nb_primes +=1
            if is_divisible_without_2(E,P,p) == False:
                n += 1
    return n/nb_primes


def density_LT_indivisible(E,P,x): # input: an elliptic curve E, a point P of E of infinite order and a real number x>0. Output: proportion of prime numbers p<x such that the Indivisibility condition holds
    N = E.conductor()
    n = 0
    nb_primes = 0
    for p in range(2,x):
        if is_prime(p) and N%p != 0:
            nb_primes +=1
            if is_divisible(E,P,p) == False:
                n += 1
    return n/nb_primes

def local_density_LT_indivisible(E,P,l,x): # input: an elliptic curve E, a point P of E of infinite order, a real number x>0 and a prime l. Output: proportion of primes p<x such that the Indivisibility condition at l holds
    n = 0
    N = E.conductor()
    nb_primes = 0
    for p in range(2,x):
        if is_prime(p) and N%p !=0:
            nb_primes += 1
            Ep = E.reduction(p)
            Pp = P.reduction(p)
            if (Ep.order() % l !=0) or Pp.division_points(l) == []:
                n += 1
    return n/nb_primes


def local_density_exponent(E,P,l,x): # input: an elliptic curve E, a point P of E of infinite order, a real number x>0 and a prime l. Output: proportion of primes p<x such that the Exponent condition at l holds
    N = E.conductor()
    n = 0
    nb_primes = 0
    for p in range(2,x):
        if is_prime(p) and (E.has_bad_reduction(p) ==False):
            nb_primes += 1
            if l_adic_valuation(((E.reduction(p)).abelian_group()).exponent(),l) == l_adic_valuation((P.reduction(p)).order(),l) :
                n += 1
    return n/nb_primes


def density_exponent_without_2(E,P,x): # input: an elliptic curve E, a point P of E of infinite order, a real number x>0 and a prime l. Output: proportion of primes p<x such that the Exponent condition holds at all primes ell different from 2
    N = E.conductor()
    n = 0
    nb_primes = 0
    for p in range(2,x):
        if is_prime(p) and N%p != 0:
            nb_primes += 1
            exponent = ((E.reduction(p)).abelian_group()).exponent()
            a = True
            for l in range(3,exponent+1):
                if is_prime(l):
                    if l_adic_valuation(exponent,l) != l_adic_valuation((P.reduction(p)).order(),l):
                        a = False
            if a:
                n += 1
    return n/nb_primes

def density_exponent(E,P,x): # input: an elliptic curve E, a point P of E of infinite order, a real number x>0 and a prime l. Output: proportion of primes p<x such that the Exponent condition holds
    N = E.conductor()
    n = 0
    nb_primes = 0
    for p in range(2,x):
        if is_prime(p) and N%p !=0:
            nb_primes += 1
            exponent = ((E.reduction(p)).abelian_group()).exponent()
            if exponent == (P.reduction(p)).order():
                n += 1
    return n/nb_primes


## CM EXAMPLE with the curve E:y^2 = x^3 - 2x and point P=E(2,2)

K.<i> = NumberField(x^2+1)
EQ = EllipticCurve([0,0,0,-2,0])
PQ = EQ(2,2)
E1 = EllipticCurve(K,[0,0,0,-2,0])
P1 = E1(2,2)



def density_LT_indivisible_CM(x): #input: a positive number x. Output: proportion of primes over the rational odd primes p<x such that the Indivisibility condition holds for the elliptic curve defined over Q(i) by the Weierstrass equation y^2 = x^3 - 2x
    n = 0
    nb_primes = 0
    for p in range(3,x):
        if is_prime(p):
            if p%4 == 3:
                nb_primes += 1
                if is_divisible_without_2(E1,P1,p) == False:
                    n += 1
            if p%4 == 1:
                nb_primes += 2
                if is_divisible_without_2(EQ,PQ,p) == False:
                    n += 2
    return n/nb_primes



## EXAMPLE with the curve E:y^2 = x^3 - 9x- 12 and P=E(-1,2)

R6 = Integers(6)
R2 = Integers(2)
R3 = Integers(3)


def rank(M,l): # returns True if the rank of the torsion part of M mod l is 0 or 1, False otherwise
    Rl = Integers(l)
    Ml = matrix(M).change_ring(Rl)
    b = False
    for k in Rl:
        if ((Ml[0][0]-1) == k * Ml[0][1] and Ml[1][0] == k * (Ml[1][1]-1) ) or (Ml[0][1] == k* (Ml[0][0]-1) and (Ml[1][1]-1) == k * Ml[1][0]):
            b = True
            return b
    return b

def last_column_in_image(M,l): # returns True if the last column of M mod l, removing the last entry, is in the image of the two first columns of (M-Id) mod l, False otherwise
    Rl = Integers(l)
    Ml = matrix(M).change_ring(Rl)
    b = False
    for k1 in Rl:
        for k2 in Rl:
            if (Ml[0][2] == k1 * (Ml[0][0]-1) + k2 * Ml[0][1] and Ml[1][2] == k1 * Ml[1][0] + k2 * (Ml[1][1] -1)):
                b = True
                return b
    return b


M1 = matrix(R6,[[1,0,0],[0,1,4],[0,0,1]])
M2 = matrix(R6,[[1,0,2],[0,1,4],[0,0,1]])
M3 = matrix(R6,[[5,0,0],[0,5,0],[0,0,1]])
M4 = matrix(R6,[[5,0,0],[5,1,0],[0,0,1]])
M5 = matrix(R6,[[1,3,0],[1,4,0],[0,0,1]])
M6 = matrix(R6,[[3,4,3],[2,3,0],[0,0,1]])
M7 = matrix(R6,[[1,4,0],[4,5,3],[0,0,1]])

G6 = MatrixGroup(M1,M2,M3,M4,M5,M6,M7) # generates the image of the mod 6 torsion-Kummer representation attached to E and P


k = 0

for g6 in G6: # counting suitable matrices in the image of the mod 6 torsion-Kummer representation for the Indivisibility condition
    if rank(g6,2) == False or last_column_in_image(g6,2) == False:
        if rank(g6,3) == False or last_column_in_image(g6,3) == False:
            k += 1

## LOCAL DENSITIES

def ind_GL2(l): # local density for the indivisibility condition if the mod l torsion-Kummer representation is surjective
    return 1 - (l^4 - 2 *l^2 - l +1)/(l^3 * (l-1) * (l^2-1))

def exp_GL2(l): # local density for the exponent condition if the ell-adic torsion-Kummer representation is surjective
    return 1- (l^5 - l^3 - l^2 -1)/(l^7 - l^6 - l^3 + l^2)


