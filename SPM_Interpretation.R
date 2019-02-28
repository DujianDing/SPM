#-------NOTE--------
# p_j,k is primarily indexed by j.
# ntree == T && nptwy == K
#-------------------

# user defined parameters
#--------------------------
nptwy = 3
ntree = 2
Tol = 0.2
in_filename = "vaf.txt"
res_filename = paste0("SPM_anc_out_K", nptwy, "_T", ntree, "_Tol", Tol, "_", in_filename)
addrin = "../DataPreprocessing/data/toy/"
addrres = "../DataPreprocessing/data/SPM_out/"
#--------------------------

vaf = read.table(paste0(addrin, in_filename), sep = '\t')
res = read.table(paste0(addrres, res_filename), sep = '\t')
res = as.matrix(res)
m = dim(vaf)[1]
n = dim(vaf)[2]

base_s = n*nptwy +m*n*nptwy
base_a = n*nptwy
base_A = n*nptwy +m*n*nptwy + m*ntree
base_p = 0

subs=rep(0,m)
pathways=rep(0,n)
for (r in 1:dim(res)[1]){
  print(paste("solution", r, "is:"))
  
  obj = 0
  for (i in 1:m) {
    for (j in 1:n) {
      for (k in 1:nptwy)
        obj = obj - res[r, base_a+(i-1)*n*nptwy+(j-1)*nptwy+k]*vaf[i,j]
    }
  }
  print(paste("objective function equals", obj))
  
  for (k in 1:nptwy){
    part = c()
    for (j in 1:n){
      if (res[r,base_p+ (j-1)*nptwy + k] == 1)
        part = c(part, j)
    }
    print(paste(c("pathway", k, "consists of gene", part), collapse = " "))
    pathways[part]=k
  }

  for (t in 1:ntree){
    part = c()
    for (i in 1:m){
      if (res[r, base_s+ (i-1)*ntree + t] == 1)
        part = c(part, i)
    }
    print(paste(c("subtype", t, "consists of", length(part), "samples:", part), collapse = " "))
    subs[part]=t
    print(paste(c("subtype", t, "is in the shape"), collapse = " "))
    for (k in 1:nptwy) {
      str = base_A+ (t-1)*nptwy*nptwy + (k-1)*nptwy + 1
      end = base_A+ (t-1)*nptwy*nptwy + k*nptwy
      print(unname(res[r,str:end]))
    }
  }
  print("-------------------")
}
