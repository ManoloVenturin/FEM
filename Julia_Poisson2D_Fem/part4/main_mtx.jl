# main_mtx.jl

using SparseArrays

# Sparse matrix - COO

row = [1, 2, 2, 3, 3]       # 1-based indexing
col = [1, 2, 4, 1, 3]
val = [10.0, 20.0, 30.0, 40.0, 50.0]

A = sparse(row, col, val)   # Automatic conversion to CSC


# Sparse matrix - CSC

rowval = [1, 3, 2, 3, 2]        # riga degli elementi (1-based)
nzval  = [10.0, 40.0, 20.0, 50.0, 30.0]  # valori non nulli
colptr = [1,    3,    4,    5,    6]     # inizio di ogni colonna
#          ↑     ↑     ↑     ↑
#        col 1 col 2 col 3 col 4

println("Sparse matrix CSC format:")
display(A)

println("rowval: ", A.rowval)
println("colptr: ", A.colptr)
println("nzval:  ", A.nzval)
