def density_LT(E,P,x): # input: an elliptic curve E, a point P of E of infinite order and a real number x>0. Output: proportion of prime numbers p<x such that P generates E(F_p)
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

def is_divisible(E,P,p): # returns True if there is a prime ell such that ell divides #Ep and P mod p is ell divisible, else returns False
    Ep = E.reduction(p)
    Pp = P.reduction(p)
    n = Ep.order()
    if is_prime(n):
        return False
    for l in range(2,n+1):
        if is_prime(l) and n%l ==0:
            if Pp.division_points(l) != []:
                return True
    return False


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


def l_adic_valuation(n,l): # input: an integer n and a prime l. Output: the l-adic valuation of n
    e = 1
    while n % (l^e) == 0:
        e += 1
    return e-1

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


def density_exponent_without_2(E,P,x): # input: an elliptic curve E, a point P of E of infinite order, a real number x>0 and a prime l. Output: proportion of primes p<x such that the Exponent condition holds at all primes l different from 2
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


## CM EXAMPLE with the curve E:y^2 = x^3 - 2x

K.<i> = NumberField(x^2+1)
EQ = EllipticCurve([0,0,0,-2,0])
PQ = EQ(2,2)
E1 = EllipticCurve(K,[0,0,0,-2,0])
P1 = E1(2,2)



def density_LT_divisible_CM(x): #input: a positive number x. Output: proportion of primes of norm <x such that the Indivisibility condition holds for the elliptic curve defined over Q(i) by the Weierstrass equation y^2 = x^3 - 2x
    n = 0
    nb_primes = 0
    for p in range(3,x):
        if is_prime(p):
            if p%4 == 3:
                nb_primes += 1
                if is_divisible_CM(E1,P1,p) == False:
                    n += 1
            if p%4 == 1:
                nb_primes += 2
                if is_divisible_CM(EQ,PQ,p) == False:
                    n += 2
    return n/nb_primes




