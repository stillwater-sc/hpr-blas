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


