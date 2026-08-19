Note taking


Okay so, there is some issue with the RBFs at the boundaries.
I will have to revist this and fix it. But before that I am going to work on, what I believe, is a significant optimization for the rbf coefficients computation.

Currently we are solving b = Ax iteratively to find x, where x are the coefficients of the rbf. 
Alternatively, if I find dx/db, that will be much more useful!

I can first put $b_{1}$ as 1 and all $b_{i\neq1}$ as 0. I can then find x for this. Then I can put $b_2$ as 0 and find x again. Since it's a linear system, when I get a new b next, it will be a linear combination of all pevious values.

This will require a lot of initial computation, especially since we'd have to do it to a very high tolerance. But it _should_ be significantly cheaper in the long run, especially if running for a long time.

I am going to try implement that for now. 


