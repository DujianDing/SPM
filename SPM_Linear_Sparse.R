library(Rcplex)
library(slam)
# user-defined parameters
C = 5

# default parameters
in_filename = "LUAD_Mutex_Matrix.txt"
K = 5

# construct constraint matrix
# Amat = read_stm_CLUTO(paste0("../DataPreprocessing/results/SPMInput/stm_K", K, "_", in_filename))
Amat = read_stm_CLUTO(paste0("stm_K", K, "_", in_filename))
order_j = order(Amat$j)
Amat$i = Amat$i[order_j]
Amat$v = Amat$v[order_j]
Amat$j = Amat$j[order_j]

#define other parameters
# M = read.table(paste0("../DataPreprocessing/data/", in_filename), sep = '\t')
M = read.table(in_filename, sep = "\t")
m = dim(M)[1]
n = dim(M)[2]
Cay = (K+1)^(K-1)

# construct the linear coefficients of objective function
obj_p = rep.int(0, times = n*K)
obj_a = rep.int(-1, times = m*K)
obj_qs = rep.int(0, times = Cay+m*Cay)

cvec = c(obj_p, obj_a, obj_qs)

# define parent indicator function h(t,k)
parind = read.table(paste0("./phylogenySet/ParentVectors_K-",K,".txt"))
h = read.table(paste0("./phylogenySet/ParentVectors_K-",K,".txt"))
hrow = dim(h)[1]
hcol = dim(h)[2]
counter = K*Cay
for (i in 1:hrow)
  for (j in 1:hcol)
    if (h[i, j] == 0) {
      h[i, j] = j
      counter = counter - 1
    }

# construct the right-hand side constraints and inequality directions
b_n_m = rep.int(1, times = n+m)
b_1 = c(C)
b_c = rep.int(0, times = Cay)
b_k = rep.int(1, times = K)
b_mk = rep.int(0, times = m*K)
b_mkc = rep.int(-1, times = m*counter)
sense_n_m = rep("E", n+m)
sense_1_c = rep("L", 1+Cay)
sense_k_mk_mkc = rep("G", K+m*K+m*counter)
bvec = c(b_n_m, b_1, b_c, b_k, b_mk, b_mkc)
sense = c(sense_n_m, sense_1_c, sense_k_mk_mkc)

# configure other parameters
vtype = "B"
objsense = "min"
maxsol = NA

# solve it by CPLEX
res = Rcplex(cvec, Amat, bvec, objsense = objsense, sense = sense, vtype = vtype, n = maxsol)

filename = 'SPM_p1_output.txt'
if (file.exists(filename))
  file.remove(filename)
lapply(res, function(x) write.table( t(x$xopt), filename  , append= T, sep='\t', row.names = F, col.names = F))

Rcplex.close()

