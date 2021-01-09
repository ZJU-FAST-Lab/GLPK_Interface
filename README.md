# GLPK_Interface
A Simple Interface for GLPK, a Linear Programming Solver for general dimensions.

# Description

linprog:
        min cTx s.t. Ax<=b
input:
        c: d*1 objective coeffs
        A: m*d constraint matrix
        b: m*1 constraint bound
        ipm: use interior point method
             or simplex method
        verbose: show details
output:
        x: d*1 decision variables
return:
        inf: No feasible solution or fail
       -inf: Unbounded problem
        real: minimum objective function

# Misc
GLPK is recommended when d the dimension of decision variable has median or 
large size. For small fixed dimension, algs that enjoy linear complexity O(m) 
are better.
