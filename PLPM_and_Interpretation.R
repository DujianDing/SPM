library(Rcplex)

# test setting
K = 2
in_filename = "toy.txt"
Amat = read_stm_CLUTO(paste0("/home/sahandk/SPM/Updates/PLPM/PLPM_K", K, "_", in_filename))
vaf = read.table(paste0("../DataPreprocessing/data/toy/", in_filename), sep = '\t')
control = list(solnpoolintensity=0, solnpoolgap=1, solnpoolagap=0, round=1)

# construct constraint matrix
order_j = order(Amat$j)
Amat$i = Amat$i[order_j]
Amat$v = Amat$v[order_j]
Amat$j = Amat$j[order_j]

M = as.data.frame((vaf>0) +0)
m = dim(M)[1]
n = dim(M)[2]

# construct the linear coefficients of objective function
obj_pjk = c()
obj_aik = rep.int(-1, times = m*K)
obj_fik = rep.int(2, times = m*K)
for (j in 1:n)
  obj_pjk = c(obj_pjk, rep(sum(M[,j]), K) )

cvec = c(obj_pjk, obj_aik, obj_fik)

# construct the right-hand side constraints and inequality directions
b_nk = rep.int(1, times = n+K)
b_mk = rep.int(0, times = m*(2*K - 1))
sense_n = rep("E", n)
sense_mk = rep("G", 2*m*K-m+K)
bvec = c(b_nk, b_mk)
sense = c(sense_n, sense_mk)

# configure other parameters
vtype = "B"
objsense = "min"
maxsol = NA

# solve it by CPLEX

res = Rcplex(cvec, Amat, bvec, objsense = objsense, sense = sense, vtype = vtype, n = maxsol, control = control)

for (r in 1:length(res)){
  print(paste("solution", r, "of cost", res[[r]]$obj, "is:"))
  for (k in 1:K){
    part = c()
    for (j in 1:n){
      if (res[[r]]$xopt[k+(j-1)*K] == 1)
        part = c(part, j)
    }
    print(paste(c("pathway", i, "consists of gene", part), collapse = " "))
  }
  print("--------")
}

# filename = paste0("PLPM_out_K", K, "_", in_filename)
# if (file.exists(filename))
#   file.remove(filename)
# invisible(lapply(res, function(x) write.table( t(x$xopt), filename  , append= T, sep='\t', row.names = F, col.names = F)))

Rcplex.close()





