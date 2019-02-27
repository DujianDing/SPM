library(Rcplex)

K = 2
in_filename = "toy.txt"
vaf = read.table(in_filename, sep = '\t')
control = list(solnpoolintensity=0, solnpoolgap=1, solnpoolagap=0, tilim=30, round=1)


M = as.data.frame((vaf>0) +0)
m = dim(M)[1]
n = dim(M)[2]

# construct the linear coefficients of objective function
sigM = integer(n)
obj_aik = rep.int(-1, times = m*K)
obj_fik = rep.int(2, times = m*K)
for (i in 1:n)
  sigM[i] = sum(M[,i])

cvec = c(rep(sigM, K), obj_aik, obj_fik)

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

# construct constraint matrix
Amat = matrix(0L, nrow = 2*m*K+K-m+n, ncol = (2*m+n)*K)

# constraint "each column is assigned to exactly one set"
basei = 0
basej = 0
for (i in 1:n){
  for (j in 1:K)
    Amat[basei+ i, basej+ (j-1)*n+i] = 1
}
# constraint "for each set Pk, at least one column is assigned to it"
basei = n
basej = 0
for (i in 1:K){
  for (j in 1:n)
    Amat[basei+ i, basej+ (i-1)*n+j] = 1
}
# constraint "for each sample the progression model is satisfied"
basei = n+K
basej = n*K
for (i in 1:m){
  for (j in 1:(K-1)){
    Amat[basei+ j + (i-1)*(K-1), basej+ j + (i-1)*K] = 1
    Amat[basei+ j + (i-1)*(K-1), basej+ j+1 + (i-1)*K] = -1
  }
}
# constraint "for each row ri, the set Pk is considered mutated if it has a 1 in ri 
# or if one of its entries in row ri is flipped to make it mutated"
basei = n+K+m*(K-1)
for (i in 1:m){
  for (j in 1:K){
    basej = 0
    for (l in 1:n){
      Amat[basei+ j + (i-1)*K, basej+ l + (j-1)*n] = M[i, l]
    }
    
    basej = n*K
    Amat[basei+ j + (i-1)*K, basej+ j + (i-1)*K] = -1
    
    basej = n*K + m*K
    Amat[basei+ j + (i-1)*K, basej+ j + (i-1)*K] = 1
  }
}

# solve it by CPLEX

res = Rcplex(cvec, Amat, bvec, objsense = objsense, sense = sense, vtype = vtype, n = maxsol, control = control)

# for (r in 1:length(res)){
#   print(paste("solution", r, "of cost", res[[r]]$obj, "is:"))
#   for (i in 1:K){
#     part = c()
#     for (j in 1:n){
#       if (res[[r]]$xopt[j+(i-1)*n] == 1)
#         part = c(part, j)
#     }
#     print(paste(c("pathway", i, "consists of gene", part), collapse = " "))
#   }
#   print("--------")
# }

filename = paste0("PLPM_out_K", K, "_", in_filename)
if (file.exists(filename))
  file.remove(filename)
invisible(lapply(res, function(x) write.table( t(x$xopt), filename  , append= T, sep='\t', row.names = F, col.names = F)))

Rcplex.close()





