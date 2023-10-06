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
from AMGXnFactorization import *

#CHANGE HERE if you want to enter the data interactively:
interactive_input = False

def computation_time(new, Xna, RFF, n, fa):
    wst1 = time.time()
    cst1 = time.process_time()
    if new:
        if Xna:
            factorization = factorization_Xn_a(RFF, n, fa)
        else:
            factorization = factorization_fXn(fa, n)
    else:
        if Xna:
            factorization = (RFF.x^n-fa).factor()
        else:
            factorization = ((fa.get_polynomial())(RFF.x^n)).factor()
    wet1 = (time.time()-wst1)  # in seconds
    cet1 = time.process_time()-cst1
    return (cet1, wet1, factorization)

def computation_time_comparison(Xna, RFF, n, fa, printing):
    (ctn, wtn, fac) = computation_time(True, Xna, RFF, n, fa)
    if printing:
        print("\t"+str(fac)+"\n\nNew AMG Algorithm: \n\tCPU time: \t"+str(ctn)+" seconds\n\tWall time: \t"+str(wtn)+" seconds\n\nSageMath Algorithm (PARI): ")

    (cto, wto, faco) = computation_time(False, Xna, RFF, n, fa)

    if not fac == faco:
        raise SystemError("The factorizations of the new and the SageMath Algorithm are not equal. Something must have gone wrong!")

    rno = int(ctn/cto)
    ron = int(cto/ctn)

    if printing:
        print("\tCPU time: \t"+str(cto)+" seconds\n\tWall time: \t"+str(wto)+" seconds\n\nratio AMG:SM = "+str(rno)+":1\nratio SM:AMG = "+str(ron)+":1")
    
    return (ctn, wtn, cto, wto, rno, ron, fac)

def str_comparison(Xna, RFF, n, fa):
    comparison = computation_time_comparison(Xna, RFF, n, fa, False)
    if Xna:
        return "\n"+str(RFF.q)+sep+str(n)+sep+str(factor(n))+sep+str(fa)+sep+str(fa.multiplicative_order())+sep+str((((comparison[6])[-1])[0]).degree())+sep+str(comparison[0])+sep+str(comparison[2])+sep+str(comparison[4])+sep+str(comparison[5])
    else:
        return "\n"+str(RFF.q)+sep+str(n)+sep+str(factor(n))+sep+str(fa.list)+sep+str(fa.get_order())+sep+str((((comparison[6])[-1])[0]).degree())+sep+str(comparison[0])+sep+str(comparison[2])+sep+str(comparison[4])+sep+str(comparison[5])

def measurements():
    name = "phdtable"
    Xnas = [(31, 675, "1"),(31, 675, 2),(31, 6075, "2"),(8, 7^4, "1"), (8, 7^4, "b"), (16, 675, "1"),(7, 648, 2), (7, 7776, 2), (7, 23328, 2),(8, 3*7^4, "1")]
    fXns = [(4, 3^4, [1,0,"b",1])] 

    filepath = "/home/uni/Unibox Rostock/Dissertation/2023-04 Entwurf/examples/"
    sep = "\t"

    table_header_Xna = "\n\nq"+sep+"n"+sep+"factor(n)"+sep+"a"+sep+"ord(a)"+sep+"max deg"+sep+"New AMG Alg"+sep+"SageMath Alg"+sep+"AMG/SM"+sep+"SM/AMG"
    table_header_fXn = "\n\nq"+sep+"n"+sep+"fator(n)"+sep+"f"+sep+"ord(f)"+sep+"max deg"+sep+"New AMG Alg"+sep+"SageMath Alg"+sep+"AMG/SM"+sep+"SM/AMG"


    file = open(filepath+"AMGXnAlg_measurements_"+name+".csv", "w")
    file.write("Computation time measurements \nfor the new X^n-a factorization algorithm \nby Anna-Maurin Graner\n\nX^n-a computations:"+table_header_Xna)

    for tupl in Xnas:
        q= tupl[0]
        RFF = RichFiniteField(q,"b","X")
        b = RFF.gen
        n = tupl[1]
        alpha = RFF.F(tupl[2])
        
        file.write(str_comparison(True,RFF,n,alpha))
    
    file.write("\n\nf(X^n) computations:"+table_header_fXn)

    for tupl in fXns:
        q= tupl[0]
        RFF = RichFiniteField(q,"b","X")
        b = RFF.gen
        n = tupl[1]
        f = RichPolynomial([RFF.F(el) for el in tupl[2]], RFF)

        file.write(str_comparison(False, RFF, n, f))

    file.close()
    print("measurements are done and can be found in "+filepath+".")
    return
    

def interactive_program():
    q=0
    while(is_prime_power(q)==False):
        q = input("\nPlease enter a prime power q (= field size): ")
        try:
            q=int(q)
        except ValueError:
            q=0

    RFF = RichFiniteField(q, "b", "X")

    print("\nThank you, we consider the following finite field: \n\t"+str(RFF))

    option = 0
    while not (option in  {1,2,3}) : 
        option = input("\nDo you wish to compute \n\t(1) X^n-a  or \n\t(2) f(X^n) for an irreducible polynomial f over Fq or \n\t(3) both? \n(1/2/3): ")
        try:
            option = int(option)
        except ValueError:
            option = 0

    n=0
    while n<1:
        n = input("\nPlease enter a positive integer n: ")
        try:
            n = int(n)
        except ValueError:
            n = 0
    
    alpha = False
    f = False 

    if option in {1,3}:
        while alpha == False:
            alpha = input("\nPlease enter an element a of Fq (e.g. b^2+1): ")
            try:
                alpha = RFF.F(alpha)
            except:
                alpha = False 
    
    if option in {2,3}:
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

    comparison = "a"
    while (not type(comparison) == bool):
        comparison = input("\nDo you wish to compare the computation of the new algorithm with the existing SageMath implementation? \n(y/n): ")
        if comparison == "n":
            comparison = False
        elif comparison == "y":
            comparison = True 
        else:
            pass                   
        
    return (option, comparison, q, RFF, n, alpha, f)

def print_factorization_time(computation):
    print("\t"+str(computation[2])+"\n\nCPU time: \t"+str(computation[0])+" seconds\nWall time: \t"+str(computation[1])+" seconds")
    return


def main():
    print("This program computes the factorization of X^n-a or f(X^n) over a finite field Fq \nwith the new formula from the paper <<The factorization of X^n-a and f(X^n)>> by Anna-Maurin Graner, \nwhere n is a positive integer and a an element of Fq.")

    #measurements()

    if interactive_input:
        (option, comparison, q, RFF, n, alpha, f) = interactive_program()
    else:
        # CHANGE HERE:
        q = 4
        n = 7^4*3

        # computation options:
        #   (1) X^n-a
        #   (2) f(X^n) for f irreducible
        #   (3) both
        option = 3

        alpha = 1
        f = [1,0,1,1] 

        # comparison with the SageMath algorithm (PARI)?
        comparison = False
    
        RFF = RichFiniteField(q, "b", "X")
        f = RichPolynomial([RFF.F(el) for el in f], RFF)
        alpha = RFF.F(alpha)

    print("\nInput is: \n\tq = "+str(q)+"\n\tFq: "+str(RFF)+"\n\n\tn = "+str(n)+" = "+str(factor(n))+"\n\n\toption = "+str(option)+"\n\tcomparison = "+str(comparison))

    if option in {1,3}:
        print("____________________\na = "+str(alpha)+"\nord(a) = "+str(alpha.multiplicative_order()))
        print("\nX^"+str(n)+"-("+str(alpha)+") = ")
        if comparison:
            computation_time_comparison(True, RFF, n, alpha, True)
        else:
            print_factorization_time(computation_time(True, True, RFF, n, alpha))

    if option in {2,3}:
        print("____________________\nf = "+str(f.get_full_info()))
        print("\nf(X^"+str(n)+") = ")
        if comparison:
            computation_time_comparison(False, RFF, n, f, True)
        else:
            print_factorization_time(computation_time(True, False, RFF, n, f))

    return 

    

if __name__ == "__main__":
    main()