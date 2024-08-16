import numpy as np
import scipy as sp

start, end = "ELSES_MATRIX_", "_A.mtx"
xxx = input("file type : ")

# OPEN INPUT FILE
infile = start + xxx + end
info   = sp.io.mminfo(infile)
matrix = sp.io.mmread(infile).tocsr()

# OPEN OUTPUT FILE
outfile = xxx + ".csr"
fp = open(outfile, 'w', encoding="utf-8")

# WRITE MATRIX INFORMATION
print(f"# {infile}", file=fp)
print(f"# {info}", file=fp)
row_ptr_size = matrix.indptr.size
col_ind_size = matrix.indices.size
data_size    = matrix.data.size
print(f"{row_ptr_size} {col_ind_size} {data_size}", file=fp)

# WRITE MATRIX DATA IN CSR FORMAT
i, j, k = 0, 0, 0
while True:
    ptr, ind, data = -1, -1, -1
    if i < row_ptr_size:
        ptr = matrix.indptr[i]
    if j < col_ind_size:
        ind = matrix.indices[j]
    if k < data_size:
        data = matrix.data[j]
    if ind == ptr == data == -1:
        break
    if info[4] == 'real':
        print("{:d} {:d} {:.20f} 0.0".format(ptr, ind, data.real), file=fp)
    if info[4] == 'complex':
        print("{:d} {:d} {:.20f} {:.20f}".format(ptr, ind, data.real, data.imag), file=fp)
    i,j,k = i+1,j+1,k+1

# CLOSE FILE
fp.close()
                                                                                                        
