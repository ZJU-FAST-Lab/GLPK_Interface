# GLPK_Interface

An Easy-to-Use GLPK Interface for General-Dimension Linear Programming.

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

GLPK is recommended for __median or large d__, i.e., the dimension of decision variable. For small dimension (d<=10), algs that enjoy linear complexity O(m) are better.
