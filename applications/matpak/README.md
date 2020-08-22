# Matrix Utility Package

MATPAK is a matrix utility package and is part of the HPRBLAS project.  Some functions make use of [MTL4](http://www.simunova.com).  

MATPAK contains replicates several standard MATLAB functions such as `fliplr`, `flipup`, `toeplitz`, `hankel`, and  `rot90`.   Other functions included are to facilitate research on centrosymmetric matrices.   

# Contents in Brief

ei - Standard basis element.  `ei(i,n)`  
fliplr - Flip matrix elements left/right  
flipud - Flip matrix up/down  
hankel - Generate hankel matrix  
isequal - determines if two matrices are equal (logical)  
ism - determines if a matrix is an M-matrix (logical)  
rot90 - Rotates a matrix 90 degrees counterclockwise  
rowsto - Generates a row stochastic matrix  
toeplitz - Generates a Toeplitz matrix    



# Contents in Detail
```
// Generate 3rd standard basis element from \R^7.
Vector ei(3,7); // Returns [0, 0, 1, 0, 0, 0, 0]
```

```
adsf
```



main()


qr(A, Matrix& Q, Matrix& R){


}
