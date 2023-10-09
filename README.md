# PolynomialsOverFiniteFields

This repository is licensed under the terms of the GNU General Public License v3.0 (GPL-3.0).

This SageMath/Python repository contains the two Python classes RichFiniteField and RichPolynomial which are useful for working with univariate polynomials over different finite fields and their extensions.

Additionally, you can find programs implementing my mathematical research which make use of these two classes.

### Table of contents:
- [Packages needed](https://github.com/amg-code/PolynomialsOverFiniteFields#packages-needed)
- [The new X^n factorization algorithm](https://github.com/amg-code/PolynomialsOverFiniteFields#the-new-xn-factorization-algorithm)
- [How rich are the two RichClasses?](https://github.com/amg-code/PolynomialsOverFiniteFields#how-rich-are-the-two-richclasses)
- [Dependencies](https://github.com/amg-code/PolynomialsOverFiniteFields#dependencies)
- [Author](https://github.com/amg-code/PolynomialsOverFiniteFields#author)

## Packages needed

The class RichPolynomial needs the Python package __multiset__. Please install the package on your system (while SageMath is running) with the command 
>sage: __pip install multiset__

You might need to restart SageMath and/or your computer afterwards. (Entering `exit` shuts down the SageMath console)


## The new X^n factorization algorithm

The file __AMGXnFactorization.py__ contains the new factorization algorithm for polynomials of the form X^n-a or f(X^n) as presented in  __"The factorization of X^n-a and f(X^n)" by Anna-Maurin Graner__ [https://arxiv.org/abs/2306.11183].

For this, 
- `factorization_Xn_a(RFF, n, alpha)`, where `RFF` is a RichFiniteField and alpha an element of the FiniteField `RFF.F`.
  This function is an exact implementation of the new factorization formula  in Theorem 18 from the paper mentioned above.
- `factorization_fXn(f,n)`, where `f` is a RichPolynomial over a RichFiniteField Fq. This function calls `factorization_Xn_a` for a root alpha of the polynomial f and the splitting field of f. Then Theorem 3 (due to Mullin,Mullen,Yucas 2010) is applied and the q-spin of all irreducible factors computed. 

The file __TestNewAlgorithm.sage__ is meant for using or testing the new factorization functions from AMGXnFactorization.py. Please search for `CHANGE HERE` in the code. These words mark the changes that the user can make to the program. The user can enter the data interactively by changing the variable `interactive_input` to `True`. Otherwise just change the parameters directly in the `main()`-function. If the variable `comparison` is set to `True`, the program first computes the factorization with the new algorithm and then with the existing SageMath-algorithm `factor()` for elements of PolynomialRings. This algorithm is based on PARI. The computation time of both algorithms and their ratio are given afterwards. 

There exists another function called `measurements()` which can be used for measuring and comparing the computation times of the new algorithm and the SageMath algorithm (PARI) for many examples at the same time without looking at the factorizations themselves. For this, please specify the `filepath` and `name` that you would like to give to your file. 

The new algorithm performs much better than the SageMath algorithm. For many integers n that are "too large" (for f(X^n) even n=81 can be too large), the SageMath algorithm either takes ages (does not return a result after a reasonable amount of time), causes the PARI stack to overflow or SageMath to crash completely. 

Some CPU computation time comparisons between the two algorithms (2023-10-09):
| q | n  | ord(a) | SageMath | AMG-Alg | ratio|
|--- | --- | --- | --- | --- | --- |
| s=1| 
| 31 | 675 = 3^3*5^2 | 1 | 0.0515 s | 0.0112 s| 4 : 1|
| 31 | 675 = 3^3*5^2 | 5 | 0.1440 s | 0.0022 s| 65 : 1|
| 31 | 6075 = 3^5*5^2 | 5 | 11.1543 s | 0.0149 s | 748 : 1| 
|31| 759,375 = 3^5*5^5 | 5 | $\infty$ s | 3.3683 s | $\infty$ : 1 |
|31 | 11,390,625 = 3^6*5^6 | 5 | $\infty$ s | 50.7298 s | $\infty$ : 1 |
| s=2 |
| 11 | 400 = 2^2*5^2 | 10|0.0570 s| 0.0093 s| 6 : 1|
| 11 | 10,000= 2^4*5^4 | 10 | 220.2496 s | 0.1791 s | 1229 : 1| 
|11| 100,000= 2^5*5^5 | 10| $\infty$ s | 3.3107 s | $\infty$ : 1 |
|11| 1,000,000 = 2^6*5^6 | 10 | Killed | 348.5296 s | $\infty$ : 1 |

Note that s is the degree of the extension field over Fq where the computations are carried out.


## How rich are the two RichClasses?
__RichFiniteField__ 

This class stores a SageMath FiniteField together with its PolynomialRing so that these two can be treated as a unit and used together. The finite field will always be of a primitive modulus so that primitive roots of unity in this finite field can easily be constructed by taking the generator to the respective exponent. 

The class has a kid called __RichExtensionField__. It stores a RichFiniteField which is an extension field of another RichFiniteField-instance - called the basefield. This class can cast elements and RichPolynomials from one of the two fields to the other. Furthermore, it can compute the minimal polynomial and the characteristic polynomial of elements of the extension field over the basefield. 

__RichPolynomial__ 

This class stores a polynomial over a given finite field (stored as a RichFiniteField) as a list and as a polynomial. This makes working with the polynomial as a list and as a polynomial in parallel easy.

It is enRICHed with many functions, some of which are redirections to existing SageMath-functions, many others are new implementations. 
All functions are split into private functions doing the computations and storing the result in a class attribute and a public function returning this attribute. This has the big advantage that all computations are done exactly once and only if needed. 

## Dependencies

The package __sage.coding.relative_finite_field_extension__ is used for the RichExtensionField and the computation of the order and the k-normality of a RichPolynomial. This package is marked as experimental. When used for the first time, it raises the following warning:
> FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
 
However, at the moment (2023-10) it only makes use of the function `Hom(F,E)` for two finite fields F and E, which is part of the standard library of SageMath.


#### Author
[Anna-Maurin Graner](https://www.mathematik.uni-rostock.de/en/struktur/professuren-apl-prof/diskrete-mathematik/translate-to-english-anna-maurin-graner/)
