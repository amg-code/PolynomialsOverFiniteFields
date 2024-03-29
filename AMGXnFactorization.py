

# This file was *autogenerated* from the file AMGXnFactorization.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_4 = Integer(4); _sage_const_2 = Integer(2)
from RichPolynomialClass import *
from RichFiniteFieldClass import *

import itertools

def k_subsets(s,k):
    return set(itertools.combinations(s,k))

# function splits the positive integer n into two positive integers n' and m such that gcd(n',q)=1 and rad(m)|q
# where q=p^l is a prime power
def make_gcd_one(q,n):
    p=((q.factor())[_sage_const_0 ])[_sage_const_0 ]
    m=_sage_const_1 
    while(n%p==_sage_const_0 ):
        n=int(n/p)
        m = m*p
    return (n, m)

# function checks whether a given integer is a prime power
def is_prime_power(q):
    if q==_sage_const_0 :
        return False
    else:
        fac = factor(q)
        if len(fac)==_sage_const_1 :
            return True
        else: 
            return False
    return 

# returns the radical of a positive integer n:
def rad(n):
    return  prod([fac[_sage_const_0 ] for fac in factor(n)])

# This function returns the cyclotomic coset of i mod d over Fq
def cycl_coset(q,d,i):
    return [i*q**j % d for j in range((Mod(q,d/(gcd(d,i))).multiplicative_order()))]

# This function returns a representative system for the q-cyclotomic classes modulo d
def cycl_coset_repsystem(q,d):
    if d==_sage_const_1 :
        coset_reps = {_sage_const_0 }
    elif d>_sage_const_1 :            
        coset_reps = set()
        tbd = list(range(_sage_const_0 ,d))
        while len(tbd)>_sage_const_0 :
            i = tbd[_sage_const_0 ]
            coset_reps.add(i)
            tbd.remove(i)
            cyclel = i*q % d 

            while not cyclel == i:
                tbd.remove(cyclel)
                cyclel = cyclel*q % d 
    return coset_reps

# function computes the q-spin of a given polynomial pol in the variable x
# coefficient degree is already given as cd 
def spin(q,x,pol,cd):
    poll = pol.list()
    spin = pol
    for j in range(_sage_const_1 ,cd):
        spin = spin * sum([poll[i]**(q**j)*x**i for i in range(len(poll))])
    return spin 

# function computes the q-spin of a binomial X^t-b given as (x,t,b)
# coefficient degree is already given as cd 
def spin_binomial(q,x,t,b,cd):
    spin =_sage_const_0 
    #print('t='+str(t)+', a='+str(a)+', cd='+str(cd))
    for i in range(cd+_sage_const_1 ):
        #print('i='+str(i))
        #print(k_subsets(range(cd),i))
        bproduct = sum([prod([b**(q**j) for j in J]) for J in k_subsets(range(cd),cd-i)])
        #print(bproduct)
        spin = spin+x**(t*i)*(-_sage_const_1 )**(cd-i)*bproduct
    #print(spin)
    return spin

# function casts polynomial in X^d to the basefield of REF
def list_to_basefield(pol, REF, x, d):
    pol = pol.list()
    pol = [pol[l] for l in range(len(pol)) if l%d ==_sage_const_0 ]
    pol = sum(REF.to_basef(pol[l])*x**(d*l) for l in range(len(pol)))
    return pol  

# function casts a polynomial over REF to the basefield of REF
def dict_to_basefield(pol, REF, x):
    pol = pol.list()
    pol = dict([(l,pol[l])  for l in range(len(pol)) if not pol[l]==REF.zero])
    new_pol = _sage_const_0 
    for l in pol:
        new_pol = new_pol + REF.to_basef(pol[l])*x**l
    return new_pol


# Based on Theorem 18, this function computes the factorization in Fqs and then converts the coefficients to the basefield 
# -> is called by factorization_Xn_a for the respective extension field F_{q^s}
def _factorization_Xn_a_inFqs(q, s, RFF, alpha, n,  e):

    # define the variables of Theorem 18:
    qs = RFF.q
    n1 = _sage_const_1 
    n2 = _sage_const_1 
    for fac in factor(n):
        if e%fac[_sage_const_0 ] ==_sage_const_0 :
            n1 = n1 * (fac[_sage_const_0 ])**(fac[_sage_const_1 ])
        else:
            n2 = n2 * (fac[_sage_const_0 ])**(fac[_sage_const_1 ])
    
    d1 = gcd(n1, int((qs-_sage_const_1 )/e))
    d2 = gcd(n2, qs-_sage_const_1 )
    
    zeta1 = RFF.get_unityroot(d1) 
    zeta2 = RFF.get_unityroot(d2)
    beta = alpha.nth_root(d1)

    #print("e = "+str(e))
    #print("n = "+str(n))
    #print("n1 = "+str(n1))
    #print("n2 = "+str(n2))
    #print("d1 = "+str(d1))
    #print("d2 = "+str(d2))
    #print("zeta1 = "+str(zeta1))
    #print("zeta2 = "+str(zeta2))
    #print("beta = "+str(beta)) 

    

    if s==_sage_const_1 :
        j_repsystem = list(range(d1))
        i_repsystem = list(range(d2))
        coeffdegs = dict([(i,_sage_const_1 ) for i in i_repsystem])
    else:
        #--------------------------------------------------------------------
        # Determine a representative system of the equivalence classes of the j-relation in Theorem 18
        if d1>_sage_const_1 :
            s1 = Mod(q, e*d1).multiplicative_order()
            #print("s1 = "+str(s1))

            j_repsystem = set()
            zeta1_tojs=dict([(zeta1**j, j) for j in range(d1)])
            #print("zeta1_tojs : "+str(zeta1_tojs))
            tbd = list(range(_sage_const_0 ,d1))
            while len(tbd)>_sage_const_0 :
                j = tbd[_sage_const_0 ]
                #print("j = "+str(j))
                j_repsystem.add(j)
                tbd.remove(j)
                for m in range(_sage_const_1 ,s1): # m=0 corresponds to j~j 
                    #print("m="+str(m))
                    #print(zeta1^(j*q^m)*beta^(q^m-1))
                    tocomparewith = zeta1**(j*q**m)*beta**(q**m-_sage_const_1 )
                    if tocomparewith == zeta1**j:
                        break
                    else:
                        tbd.remove(zeta1_tojs[zeta1**(j*q**m)*beta**(q**m-_sage_const_1 )])
        else:
            j_repsystem = {_sage_const_0 }

        #------------------------------------------------- 
        # representative system for selection of is:
        i_repsystem = cycl_coset_repsystem(q,d2)
        s2 = Mod(q,d2).multiplicative_order()
        #print("s2 = "+str(s2))

        tis = dict()

        ts = s2.divisors()
        ts.sort()
        d2sd2t = dict([(t,int(d2/gcd(d2, q**t-_sage_const_1 ))) for t in ts])

        for i in i_repsystem:
            for l in range(len(ts)):
                if i%d2sd2t[ts[l]]==_sage_const_0 :
                    tis[i] = ts[l]
                    break
    
    #print("j_reps : "+str(j_repsystem))
    #print("i_reps : "+str(i_repsystem))
    #print("tis : "+str(tis))
    
    
    #-----------------------------------------------------------------

    r = n2.inverse_mod(e*d1)
    n1d1 = int(n1/d1)
    irred_factors = []   
    if s>_sage_const_1 :
        X = (RFF.get_basefield()).x

    # Product over j, v, i of Theorem 18
    for j in j_repsystem:
        for v in divisors(int(n2/d2)):
            for i in i_repsystem:
                if gcd(i,v)==_sage_const_1 :
                    for m in range(gcd(s1,tis[i])): 
                        #print("(j,v,i,m)=("+str(j)+","+str(v)+","+str(i)+","+str(m)+")")                   
                        #R_jvi = RFF.x^(n1d1*v)-zeta2^i*(zeta1^j*beta)^(r*v)
                            
                        fact = spin_binomial(q,RFF.x,n1d1*v,zeta2**(i*q**m)*(zeta1**j*beta)**(r*v),lcm(s1,tis[i]))

                        if s>_sage_const_1 :
                            fact = list_to_basefield(fact, RFF, X, n1d1*v)
                        
                        irred_factors.append((fact,_sage_const_1 ))
    return irred_factors               


def factorization_Xn_a(RFF, n, alpha):
    ptol = _sage_const_1 

    # Case gcd(q,n)>1: Using remark of Chapter 2 we consider b instead of alpha where b^(p^l)=alpha
    if gcd(RFF.q, n)>_sage_const_1 :
        (n, ptol) = make_gcd_one(RFF.q, n)
        #print("n = "+str(n)+"*"+str(ptol))
        alpha = alpha.nth_root(ptol)
        #print("b = "+str(alpha)+" satisfies b^"+str(ptol)+"="+str(alpha^ptol))


    e = alpha.multiplicative_order()

    # Determine the extension field such that rad(n)|(q^s-1) and 4 !| n or q^s = 1 mod 4
    w = Mod(RFF.q,rad(n)).multiplicative_order()
    #print("rad(n)|q^w-1 for w="+str(w))
    if (not n%_sage_const_4 ==_sage_const_0 ) or RFF.q**w%_sage_const_4  ==_sage_const_1 :
        s = w
    else:
        s = _sage_const_2 *w
    #print("s = "+str(s))
    #print("q^s-1 = "+str(RFF.q^s-1)+" = "+str(factor(RFF.q^s-1)))

    if s>_sage_const_1 :
        REF = RichExtensionField(RFF, s , "y", "T")        
        irred_factors = _factorization_Xn_a_inFqs(RFF.q, s, REF, REF.to_extf(alpha), n,  e)
    else:
        irred_factors = _factorization_Xn_a_inFqs(RFF.q, s, RFF, alpha, n,  e)

    #irred_factors = list(irred_factors.items())
    if ptol>_sage_const_1 :
        irred_factors=[(fact[_sage_const_0 ], fact[_sage_const_1 ]*ptol) for fact in irred_factors]

    return Factorization(irred_factors)

def factorization_fXn(f,n):
    if not f.is_irreducible():
        raise ValueError("The given polynomial f is not irreducible.")
    else:
        REF = f.get_splitting_field()
        #print(REF)
        #alpha = (f.get_MP_in_splitting_field()).root()
        #print(alpha)
        Xnalpha = factorization_Xn_a(REF, n , (f.get_RP_in_splitting_field()).root())
        #print(Xnalpha)
        return Factorization([(dict_to_basefield(spin((f.get_RFF()).q, REF.x, fact[_sage_const_0 ], REF.extension_degree()), REF, (REF.get_basefield()).x), fact[_sage_const_1 ]) for fact in Xnalpha])
    return

