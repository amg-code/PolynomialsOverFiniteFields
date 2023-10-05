# PolynomialsOverFiniteFields
This SageMath/Python repository contains the two Python classes RichFiniteField and RichPolynomial which are useful for working with univariate polynomials over different finite fields and their extensions.

Additionally, you can find programs implementing my mathematical research which make use of these two classes.

## Packages

The class RichPolynomial needs the Python package __multiset__. Please install the package on your system (while Sage is running) with the command 
>__pip install multiset__ 

You might need to restart Sage and/or your computer afterwards.

Also the package __sage.coding.relative_finite_field_extension__ is used for the computation of the order and the k-normality of the polynomials. This package is marked as experimental. When used for the first time, it raises the following  
> FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
 
However, at the moment (2023-10) it only makes use of the function Hom(F,E) for two finite fields F and E, which is part of the standard library of SageMath.

## AMG-NewXnFactorizationAlgorithm.sage

This program is an exact implementation of the new factorization formula for polynomials of the form X^n-a over a finite field in Theorem 18 from 

__"The factorization of X^n-a and f(X^n)" by Anna-Maurin Graner__ [https://arxiv.org/abs/2306.11183].

Furthermore, it makes use of this new algorithm to derive the factorization of polynomials of the form f(X^n) where  f is an irreducible polynomial over a finite field.

Please search for _CHANGE HERE_ in the code. These words mark the changes that the user can make to the program. 

The user can either enter the data interactively by changing the variable `interactive_input` to `True`. Otherwise just change the parameters directly in the `main()`-function. 

The program first computes the factorization with the new algorithm and then with the existing SageMath-algorithm `factor()` for elements of PolynomialRings. This algorithm is based on PARI. The wall and the CPU computation time of both algorithms are given and their ratio computed. 

The new algorithm performs much better for all positive integers that are not sufficiently small. For many integers n that are "too large", the PARI stack overflows or SageMath crashes completely. 


## How rich are the two RichXXClasses?
__RichFiniteField__ 

This class stores a SageMath FiniteField together with its PolynomialRing so that these two can be treated as a unit and used together. The finite field will always be of a primitive modulus so that primitive roots of unity in this finite field can easily be constructed by taking the generator to the respective exponent. 

The class has a kid called __RichExtensionField__. It stores a RichFiniteField which is an extension field of another RichFiniteField-instance - called the basefield. This class can cast elements and RichPolynomials from one of the two fields to the other. Furthermore, it can compute the minimal polynomial and the characteristic polynomial of elements of the extension field over the basefield. 

__RichPolynomial__ 

This class stores a polynomial over a given finite field (stored as a RichFiniteField) as a list and as a polynomial. This makes working with the polynomial as a list and as a polynomial in parallel easy.

It is enRICHed with many functions, some of which are redirections to existing SageMath-functions, many others are new implementations. 
All functions are split into private functions doing the computations and storing the result in a class attribute and a public function returning this attribute. This has the big advantage that all computations are done exactly once and only if needed. 


#### Author
[Anna-Maurin Graner](https://www.mathematik.uni-rostock.de/en/struktur/professuren-apl-prof/diskrete-mathematik/translate-to-english-anna-maurin-graner/)
