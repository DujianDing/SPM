library(Rcplex)
library(slam)
#-------NOTE--------
# p_j,k is primarily indexed by j.
# ntree == T && nptwy == K
#-------------------
# user-defined parameters
in_filename = "TCGA_COAD_CCFs_spm.txt"
ntree = 1
nptwy = 3
Tol=0.2
control = list(solnpoolintensity=0, solnpoolgap=1, solnpoolagap=0, tilim=108000, round=1)

# test setting
# Amat = read_stm_CLUTO(paste0("../DataPreprocessing/results/SPMInput/ccs_K", nptwy, "_T", ntree, "_", in_filename))
# vaf = read.table(paste0("../DataPreprocessing/data/COMB/", in_filename), sep = '\t')

# production setting
Amat = read_stm_CLUTO(paste0("/home/sahandk/SPM/COAD/COAD_CCS/ccs_K", nptwy, "_T", ntree, "_", in_filename))
vaf = read.table(paste0("/home/sahandk/SPM/COAD/",in_filename), sep = "\t")

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
obj_sA = rep.int(0, times = m*ntree+ntree*(nptwy+1)*(nptwy+1))

cvec = c(obj_p, obj_alpha, obj_sA)

# construct the right-hand side constraints and inequality directions
b_punq = rep.int(1, times = n)
b_pnemp = rep.int(1, times = nptwy)
b_sunq = rep.int(1, times = m)
b_snemp = rep.int(1, times = ntree)
b_cndmut = rep.int(0, times = m*n*nptwy)
b_mutex = rep.int(1, times = m*nptwy)
b_lprg = rep.int(1, times = 2*ntree*(nptwy+1))
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
s_lprg = rep("E", 2*ntree*(nptwy+1))
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

# base_s = n*nptwy +m*n*nptwy
# base_a = n*nptwy
# base_A = n*nptwy +m*n*nptwy + m*ntree
# base_p = 0
# if (length(res) > 0){
#   
#   for (r in 1:length(res)){
#     print(paste("solution", r, "is:"))
#     print(paste("optimal cost is", res[[r]]$obj))    
#     
#     for (k in 1:nptwy){
#       part = c()
#       for (j in 1:n){
#         if (res[[r]]$xopt[base_p+ (j-1)*nptwy + k] == 1)
#           part = c(part, j)
#       }
#       print(paste(c("pathway", k, "consists of gene", part), collapse = " "))
#     }    
#     
#     for (t in 1:ntree){
#       part = c()
#       for (i in 1:m){
#         if (res[[r]]$xopt[base_s+ (i-1)*ntree + t] == 1)
#           part = c(part, i)
#       }
#       print(paste(c("subtype", t, "consists of samples", part), collapse = " "))
#       print(paste(c("subtype", t, "is in the shape"), collapse = " "))
#       for (k in 1:(nptwy+1)) {
#         str = base_A+ (t-1)*(nptwy+1)*(nptwy+1) + (k-1)*(nptwy+1) + 1
#         end = base_A+ (t-1)*(nptwy+1)*(nptwy+1) + k*(nptwy+1)
#         print(res[[r]]$xopt[str:end])
#       }
#     }    
#     print("-------------------")
#   }
# } else
#   print("No solution found!")


filename = paste0("/home/sahandk/SPM/COAD/SPM_ccs_out_K",nptwy,"_T",ntree,"_Tol",Tol,"_", in_filename)
if (file.exists(filename))
  file.remove(filename)
invisible(lapply(res, function(x) write.table( t(x$xopt), filename  , append= T, sep='\t', row.names = F, col.names = F)))


Rcplex.close()





