library(Rcplex)
library(slam)
#-------NOTE--------
# p_j,k is primarily indexed by j.
# ntree == T && nptwy == K

#-------------------
# user-defined parameters
ntree = 2
nptwy = 3
Tol=0
in_filename = "vaf.txt"
int_filename = paste0("anc_K", nptwy, "_T", ntree, "_", in_filename)
out_filename = paste0("SPM_anc_out_K",nptwy,"_T",ntree,"_Tol",Tol,"_", in_filename)
addrin = "/home/sahandk/SPM/COAD/"
addrint = "/home/sahandk/SPM/COAD/COAD_Ancestry/"
addrout = "/home/sahandk/SPM/COAD/"
control = list(solnpoolintensity=0, solnpoolgap=1, solnpoolagap=0, tilim=108000, round=1)
#-------------------


# production setting
vaf = read.table(paste0(addrin, in_filename), sep = "\t")
Amat = read_stm_CLUTO(paste0(addrint, int_filename))

# construct constraint matrix
order_j = order(Amat$j)
Amat$i = Amat$i[order_j]
Amat$v = Amat$v[order_j]
Amat$j = Amat$j[order_j]

#define other parameters
M = as.data.frame((vaf>0) +0)
m = dim(M)[1]
n = dim(M)[2]

# construct the linear coefficients of objective function
obj_p = rep.int(0, times = n*nptwy)
obj_alpha = c()
for (i in 1:m)
  for (j in 1:n)
    obj_alpha = c(obj_alpha, rep(-vaf[i,j], nptwy))
obj_sA = rep.int(0, times = m*ntree+ntree*nptwy*nptwy)

cvec = c(obj_p, obj_alpha, obj_sA)

# construct the right-hand side constraints and inequality directions
b_punq = rep.int(1, times = n)
b_pnemp = rep.int(1, times = nptwy)
b_sunq = rep.int(1, times = m)
b_snemp = rep.int(1, times = ntree)
b_cndmut = rep.int(0, times = m*n*nptwy)
b_mutex = rep.int(1, times = m*nptwy)
b_lprg = rep.int(nptwy-1, times = ntree*nptwy)
b_prg = rep.int(-2-Tol, times = m*ntree*nptwy*(nptwy-1))
b_nfng = c()
for (i in 1:m)
  for (j in 1:n)
    b_nfng = c(b_nfng, M[i,j])
s_punq = rep("E", n)
s_pnemp = rep("G", nptwy)
s_sunq = rep("E", m)
s_snemp = rep("G", ntree)
s_cndmut = rep("L", m*n*nptwy)
s_mutex = rep("L", m*nptwy)
s_lprg = rep("E", ntree*nptwy)
s_prg = rep("G", m*ntree*nptwy*(nptwy-1))
s_nfng = rep("L", m*n)

bvec = c(b_punq, b_pnemp, b_sunq, b_snemp, b_cndmut, b_mutex, b_lprg, b_prg, b_nfng)
sense = c(s_punq, s_pnemp, s_sunq, s_snemp, s_cndmut, s_mutex, s_lprg, s_prg, s_nfng)

# configure other parameters
vtype = "B"
objsense = "min"
maxsol = NA

# solve it by CPLEX
res = Rcplex(cvec, Amat, bvec, objsense = objsense, sense = sense, vtype = vtype, n = maxsol, control = control)

filename = paste0(addrout, out_filename)
if (file.exists(filename))
  file.remove(filename)
invisible(lapply(res, function(x) write.table( t(x$xopt), filename  , append= T, sep='\t', row.names = F, col.names = F)))


Rcplex.close()





