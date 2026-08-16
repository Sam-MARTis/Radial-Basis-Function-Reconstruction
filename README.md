
TODO:

1. Integration
So for 1D we can write a closed form integration for the RBFs.
For 2D and higher order, we need to do gaussian quadrature or hemite integration. 

We first get the rbfs-cell averages relations of the form:
u = Ax, where u are the vector of cell averages, x is the value of the rbs, A is the coefficients of interactions thata ccounts for the positions of the RBFs in the domain.

A is constant throughout the simulation since the RBFs' centers do not move. it would thus be worthwhile to take an inverse of the matrix rather than use an iterative method. 

Alternatively, we could find partial(x)/partial(u). Since it is linear, I just need to do 1 base computation to find x, then do N calculations of u_modified, where N is the number of cells (also the dimension of the x and u vector). u_modified (s) will be computed via just incrementing/decrementing u in one place by a small amount.





2. 

