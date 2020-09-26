

// A \in M(p,q) and B \in M(m,n), then kron(A,B) \in M(mxp, nxq)

/*
kron(X,Y) is the Kronecker tensor product of X and Y.
    The result is a large matrix formed by taking all possible
    products between the elements of X and those of Y. For
    example, if X is 2 by 3, then kron(X,Y) is
 
       [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
         X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
 
    If either X or Y is sparse, only nonzero elements are multiplied
    in the computation, and the result is sparse.
*/