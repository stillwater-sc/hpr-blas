# Matrix Utility Package

MATPAK is a collection of free functions of matrix utilities.  It is included in the HPRBLAS project.  Some functions make use of [MTL4](http://www.simunova.com).    

MATPAK contains replicates several standard MATLAB functions such as `fliplr`, `flipup`, `toeplitz`, `hankel`, and  `rot90`.   The project was originally developed to facilitate research on row stochastic centrosymmetric matrices.   

# Contents in Brief

Let A, B be matrices and v a vector.

bands - extracts diagonals from a matrix (e.g., tridiagonal)
diag - creates a diagonal matrix or extracts the diagonal
diam(A) - calculates diameter of a matrix
ek - Standard basis element.  `ek(i,n)`  
fliplr - Flip matrix elements left/right  
flipud - Flip matrix up/down
gt - greater than, e.g., A > B, returns 1/0 matrix  
hankel - Generate hankel matrix  
iscentro - determines if A is centrosymmetric (bool)
isequal - determines if two matrices are equal (bool)  
ism - determines if a matrix is an M-matrix (bool)  
isnormal - determines if a matrix A is normal
kron(A,B) - Kronecker product of A with B
mkcentro - makes a centrosymmetric matrix
reshape -  Maps m x n --> r x c where mn = rc
rot90 - Rotates a matrix 90 degrees counterclockwise  
rowsto - Generates a row stochastic matrix  
size - caclulates the size of a matrix
sum - computes sum of entries per dimension
toeplitz - Generates a Toeplitz matrix (see hankel)    



# Contents in Detail

## any

## diag
Returns the diagonal of matrix A.  If A is a vector (i.e., n x 1 matrix), diag(A) returns a n x n diagonal matrix with diagonal elements A.  

``` 
diag(A)
```


## diam
Returns the diameter of a matrix.  

``` 
diam(A)
```



## ek
Returns and elementary vector (0, 0, 0, ..., 1, 0, ..., 0)
``` 
// k = index containing 1
// n = length of vector
ek(k, n)
```
 


## Hankel
Returns a Hankel matrix.
hankel(c,r) is a non-symmetric Toeplitz matrix having C as its
    first column and R as its first row rotated 90 degrees.

```
hankel(c,r)
```

-----------------------

# ISA

## iscentro
Determines if a matrix is centrosymmetric.  

``` 
iscentro(A)
``` 

## isdiagdom
Determines if a matrix is diagonally dominate.  

``` 
isdiagdom(A)
``` 

## isequal
Determines which elements are equal in two m x n matrices, A and B. 

``` 
isequal(A,B)
``` 

## ism
Determines if a matrix is an M-matrix.  
An M-matrix is a Z-matrix (non-positive off-diagonal entries) with real part of eigenvalues nonnegative.
There are several equivalent definitions of M-matrix.  Most spatial finite difference discretizations yield M-matrices.  Moreover, M-matrix guarantees the monotonicity of the solution (i.e., non-oscillations in discrete solution, i.e., Maximum Principle preserved). See Varga,R.S. (1966). On a discrete maximum principle. SIAM J.Numer.Anal. 3, 355–359.

``` 
ism(A)
``` 

## isnormal
Determines if a matrix is *normal*: A'A = AA'.  Returns boolean.

``` 
isnormal(A)
``` 

---------------------



## rot90
Rotate array 90 degrees counterclockwise.  

``` 
rot90(A)
``` 


## rowsto
Generates a random m x n row stochastic matrix. 

``` 
rowsto(m,n)
``` 




## size
Determines if a square matrix A is centrosymmetric.  Returns boolean.
``` 
iscentro(A)
``` 



## sum 
Given a matrix A, sum(A,d) sums columns (d=1), rows (d=2), or all (d=0) entries.

```
sum(A,d)
```




## Toeplitz
Returns a Toeplitz matrix.
toeplitz(c,r) is a non-symmetric Toeplitz matrix having C as its
    first column and R as its first row.

``` 
toeplitz(c,r)
``` 

 


 
