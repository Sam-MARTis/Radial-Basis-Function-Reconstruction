Note taking


Okay so, there is some issue with the RBFs at the boundaries.
I will have to revist this and fix it. But before that I am going to work on, what I believe, is a significant optimization for the rbf coefficients computation.

Currently we are solving b = Ax iteratively to find x, where x are the coefficients of the rbf. 
Alternatively, if I find dx/db, that will be much more useful!

I can first put $b_{1}$ as 1 and all $b_{i\neq1}$ as 0. I can then find x for this. Then I can put $b_2$ as 0 and find x again. Since it's a linear system, when I get a new b next, it will be a linear combination of all pevious values.

This will require a lot of initial computation, especially since we'd have to do it to a very high tolerance. But it _should_ be significantly cheaper in the long run, especially if running for a long time.

I am going to try implement that for now. 


---

Okay so, we need a complete local interpolation instead of a global one with compact support. We can move onto the Euler equation directly. The vector will have to obtimized and a more dedicated solver instead of the current gauss siedel version would have to be used.
Might not even be needed if nearest neighbours based interpolation is used.

Let us take M nearest neighbours, with the RBFs at the cell centers. 
Then we can find the structural metrix Ai for each cell. I will then compute A inv numerically. 

Since A is a small matrix of size MxM, it should be cheap to compute the inverse. 

Then we can find the coefficients of the RBFs for each cell by doing x = A_inv * b, where b is the vector of cell averages and x are the rbf coefficients.

Ok, so. Setup:
- Loop through all cells and do:
  - Setup the M nearest neighbour indices
  - Compute the structural matrix Ai, unique to eacl cell and dependent on the local mesh structure
  - Compute the inverse of Ai
- Loop through all cells again:
  - Setup initial conditions and fluxes.

And the loop: 
- Loop through all cells again:
  - Compute and rbf coefficiencts for each coefficient of Ui by using the inverse of Ai and the cell averages.
  - Recontruct the solution to the cell edges using the rbf coefficients and the rbf basis functions.
  - Push the value to the edges. 
- Apply boundary conditions
  - Fill the missing pushed value for edges connected to just one cell. 
- Loop through all edges. 
  - Compute the numerical flux value for each edge using the pushed value.
  - Update cell averages. 


Let us skip the linear advection and move directly to the Euler equation in 1D. 

Roe flux should be sufficient


So for what class does what. 
The rbf computation should probably be handled by the mesh itself. It has all the necesary information and the solver can just query the rbf coefficients. 


The neighbourhood struct should include the current index as well in the indices of neighbourhood. 
This will be better for uniformity since it has to store the rbf coefficient for the current cell as well. 
