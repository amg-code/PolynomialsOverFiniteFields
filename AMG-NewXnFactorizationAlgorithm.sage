# This program computes the factorization of X^n-a for any positive integer n and any element a of a finite field Fq over this finite field Fq
# For this it uses the new formula from the paper <<The factorization of X^n-a and f(X^n)>> by Anna-Maurin Graner

# ****************************************************************************
#       Copyright (C) 2023 Anna-Maurin Graner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import time 
from RichPolynomialClass import *
from RichFiniteFieldClass import *

#CHANGE HERE if you want to enter the data interactively:
interactive_input = False

# function checks whether a given integer is a prime power
def is_prime_power(q):
    if q==0:
        return False
    else:
        fac = factor(q)
        if len(fac)==1:
            return True
        else: 
            return False
    return 

# returns the radical of a positive integer n:
def rad(n):
    return  prod([fac[0] for fac in factor(n)])

# This function returns a representative system for the q-cyclotomic classes modulo d
def cycl_coset_repsystem(q,d):
    if d==1:
        coset_reps = {0}
    elif d>1:            
        coset_reps = set()
        tbd = list(range(0,d))
        while len(tbd)>0:
            i = tbd[0]
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
    for j in range(1,cd):
        spin = spin * sum([poll[i]^(q^j)*x^i for i in range(len(poll))])
    return spin 

# function casts polynomial in X^d to the basefield of REF
def list_to_basefield(pol, REF, x, d):
    pol = pol.list()
    pol = [pol[l] for l in range(len(pol)) if l%d ==0]
    pol = sum(REF.to_basef(pol[l])*x^(d*l) for l in range(len(pol)))
    return pol  

# function casts a polynomial over REF to the basefield of REF
def dict_to_basefield(pol, REF, x):
    pol = pol.list()
    pol = dict([(l,pol[l])  for l in range(len(pol)) if not pol[l]==REF.zero])
    new_pol = 0
    for l in pol:
        new_pol = new_pol + REF.to_basef(pol[l])*x^l
    return new_pol


# Based on Theorem 18, this function computes the factorization in Fqs and then converts the coefficients to the basefield 
# -> is called by factorization_Xn_a for the respective extension field F_{q^s}
def _factorization_Xn_a_inFqs(q, s, RFF, alpha, n,  e):

    # define the variables of Theorem 18:
    qs = RFF.q
    n1 = 1
    n2 = 1
    for fac in factor(n):
        if e%fac[0] ==0:
            n1 = n1 * (fac[0])^(fac[1])
        else:
            n2 = n2 * (fac[0])^(fac[1])
    
    d1 = gcd(n1, int((qs-1)/e))
    d2 = gcd(n2, qs-1)
    
    zeta1 = RFF.get_unityroot(d1) 
    zeta2 = RFF.get_unityroot(d2)

    #------------------------------------------------
    # Find the element beta such that beta^d1 = alpha
    if d1>1:
        for gamma in RFF.F:
            if gamma^d1 == alpha:
                beta = gamma
                break
    elif d1==1:
        beta = alpha
    else:
        raise ValueError("Something went very wrong. Did you choose a negative value for n?")

    if s==1:
        j_repsystem = list(range(d1))
        i_repsystem = list(range(d2))
        coeffdegs = dict([(i,1) for i in i_repsystem])
    else:
        #--------------------------------------------------------------------
        # Determine a representative system of the equivalence classes of the j-relation in Theorem 18
        if d1>1:
            j_repsystem = set()
            zeta1_tojs=dict([(zeta1^j, j) for j in range(d1)])
            tbd = list(range(0,d1))
            while len(tbd)>0:
                j = tbd[0]
                j_repsystem.add(j)
                tbd.remove(j)
                for m in range(1,s): # m=0 corresponds to j~j 
                    tocomparewith = zeta1^(j*q^m)*beta^(q^m-1)
                    if tocomparewith == zeta1^j:
                        break
                    else:
                        tbd.remove(zeta1_tojs[zeta1^(j*q^m)*beta^(q^m-1)])
        else:
            j_repsystem = {0}

        #------------------------------------------------- 
        # q-cyclotomic coset mod d2 for selection of is:
        i_repsystem = cycl_coset_repsystem(q,d2)

        #---------------------------------------------------------------
        # Determine the q-coefficient degree of the irreducible factors
        s1 = Mod(q, e*d1).multiplicative_order()
        if s==s1:
            coeffdegs = dict([(i,s) for i in i_repsystem])
        else:
            coeffdegs = dict()
            # candidates for t:
            ts = [s1*l for l in range(1,int(s/s1)) if s%(s1*l)==0]
            # fractions d_2^(s)/d_2^(t) for all t in ts:
            d2d2ts = dict([(t, int(d2/gcd(n2, q^t-1))) for t in ts])
            ts.append(s)
            d2d2ts[s]=1

            # Select minimum of ts such that fraction divides i
            for i in i_repsystem:
                coeffdegs[i] = min([t for t in ts if i%d2d2ts[t]==0])  
    
    #-----------------------------------------------------------------

    r = n2.inverse_mod(e*d1)
    n1d1 = int(n1/d1)
    irred_factors = dict()    
    X = (RFF.get_basefield()).x

    # Product over j, v, i of Theorem 18
    for j in j_repsystem:
        for v in divisors(int(n2/d2)):
            for i in i_repsystem:
                if gcd(i,v)==1:                        
                    fact = spin(q, RFF.x, RFF.x^(n1d1*v)-zeta2^i*(zeta1^j*beta)^(r*v), coeffdegs[i])
                    if s>1:
                        fact = list_to_basefield(fact, RFF, X, n1d1*v)                    

                    try: 
                        irred_factors[fact] = irred_factors[fact]+1
                    except KeyError:
                        irred_factors[fact] = 1 
                    
    return Factorization(list(irred_factors.items()))


def factorization_Xn_a(RFF, n, alpha):
    e = alpha.multiplicative_order()

    # Determine the extension field such that rad(n)|(q^s-1) and 4 !| n or q^s = 1 mod 4
    w = Mod(RFF.q,rad(n)).multiplicative_order()
    #print("rad(n)|q^w-1 for w="+str(w))
    if (not n%4==0) or RFF.q^w%4 ==1:
        s = w
    else:
        s = 2*w

    if s>1:
        REF = RichExtensionField(RFF, s , "y", "T")        
        return _factorization_Xn_a_inFqs(RFF.q, s, REF, REF.to_extf(alpha), n,  e)
    else:
        return  _factorization_Xn_a_inFqs(RFF.q, s, RFF, alpha, n,  e)

    return 

def factorization_fXn(f,n):
    if not f.is_irreducible():
        raise ValueError("The given polynomial f is not irreducible.")
    else:
        REF = f.get_splitting_field()
        Xnalpha = factorization_Xn_a(REF, n , (f.get_RP_in_splitting_field()).root())
        
        return Factorization([(dict_to_basefield(spin((f.get_RFF()).q, REF.x, fact[0], REF.extension_degree()), REF, (REF.get_basefield()).x), fact[1]) for fact in Xnalpha])

def main():
    print("This program computes the factorization of X^n-a or f(X^n) over a finite field Fq \nwith the new formula from the paper <<The factorization of X^n-a and f(X^n)>> by Anna-Maurin Graner, \nwhere n is a positive integer and a an element of Fq.")

    if interactive_input:
        q=0
        while(is_prime_power(q)==False):
            q = input("\nPlease enter a prime power q (= field size): ")
            try:
                q=int(q)
            except ValueError:
                q=0

        RFF = RichFiniteField(q, "b", "X")

        print("\nThank you, we consider the following finite field: \n\t"+str(RFF))

        comp_fXn = "a"
        while (not type(comp_fXn) is bool) : 
            comp_fXn = input("\nDo you wish to compute \n\t(1) X^n-a  or \n\t(2) f(X^n) for an irreducible polynomial f over Fq? \n(1/2): ")
            if comp_fXn == "1":
                comp_fXn = False
            elif comp_fXn == "2":
                comp_fXn=True 
            else:
                pass

        if comp_fXn:
            f = False
            while f == False:
                try:
                    f  = list((((input("\nPlease enter the coefficient vector of an irreducible polynomial over Fq. \n(eg. [1,0,b,b^2] for X^3+bX+b^2)\nf = "))[1:-1]).strip()).split(","))
                    f= RichPolynomial([RFF.F(el) for el in f], RFF)
                    if not f.is_irreducible():
                        print("\nError: The polynomial "+str(f.get_polynomial())+" is not irreducible over F_{"+str(RFF.q)+"}.")
                        f = False 
                    else: 
                        print("\n f = "+str(f.get_full_info()))
                except:
                    f= False
            n=0
            while n<1:
                n = input("\nPlease enter a positive integer n: ")
                try:
                    n = int(n)
                except ValueError:
                    n = 0                
        else: 
            alpha=False
            while alpha == False:
                alpha = input("\nPlease enter an element a of Fq (e.g. b^2+1): ")
                try:
                    alpha = RFF.F(alpha)
                except:
                    alpha = False 

            n=0
            while n<1:
                n = input("\nPlease enter a positive integer n: ")
                try:
                    n = int(n)
                except ValueError:
                    n = 0
    else:
        q=8
        RFF = RichFiniteField(q,"b","X")
        b = RFF.gen
        comp_fXn = False
        f = [1,0,b,1] 
        f = RichPolynomial(f, RFF)
        alpha = RFF.F(b)
        n = 7^4*3  
        print("\nInput is: \n\tq = "+str(q)+"\n\tFq = "+str(RFF)+"\n\tn = "+str(n))

    if comp_fXn==False:
        print("\ta = "+str(alpha)+" \n \tord(a) = "+str(alpha.multiplicative_order()))
        print("\nNew algorithm by AMG: \n X^"+str(n)+"-("+str(alpha)+") =")
        wst1 = time.time()
        cst1 = time.process_time()
        f1 = factorization_Xn_a(RFF,n,alpha)
        wet1 = (time.time()-wst1)  # in seconds
        cet1 = time.process_time()-cst1
        print(str(f1)+"\n\nCPU time "+str(cet1)+" seconds\nWall time "+str(wet1)+" seconds") 

        print("\nSageMath-Algorithm: \n X^"+str(n)+"-("+str(alpha)+") =")
        wst2 = time.time()
        cst2 = time.process_time()
        f2 = (RFF.x^n-alpha).factor()
        wet2 = (time.time()-wst2)
        cet2 = time.process_time()-cst2
        print(str(f2)+"\n\nCPU time "+str(cet2)+" seconds\nWall time "+str(wet2)+" seconds")  

        print("ratio of Sagemath-Alg/AMG-Alg  (CPU): "+str(cet2/cet1) )
        print("ratio of AMG-Alg/Sagemath-Alg (CPU): "+str(cet1/cet2))   
    else: 
        print("\nf = "+str(f.get_full_info())+"\n\nNew algorithm by AMG: \n f(X^"+str(n)+") =")   
        wst1 = time.time()
        cst1 = time.process_time()
        f1 = factorization_fXn(f, n)
        wet1 = (time.time()-wst1)  # in seconds
        cet1 = time.process_time()-cst1
        print(str(f1)+"\n\nCPU time "+str(cet1)+" seconds\nWall time "+str(wet1)+" seconds")   

        print("\nSageMath-Algorithm: \n f(X^"+str(n)+") =")
        wst2 = time.time()
        cst2 = time.process_time()
        f2 = ((f.get_polynomial())(RFF.x^n)).factor()
        wet2 = (time.time()-wst2)
        cet2 = time.process_time()-cst2
        print(str(f2)+"\n\nCPU time "+str(cet2)+" seconds\nWall time "+str(wet2)+" seconds")

        print("ratio of Sagemath-Alg/AMG-Alg  (CPU): "+str(cet2/cet1) )
        print("ratio of AMG-Alg/Sagemath-Alg (CPU): "+str(cet1/cet2))  
    return

if __name__ == "__main__":
    main()