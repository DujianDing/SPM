#-------NOTE--------
# p_j,k is primarily indexed by j.
# ntree == T && nptwy == K
#-------------------

# user defined parameters
#--------------------------
in_filename = "TCGA_COAD_CCFs_spm.txt"
nptwy = 5
ntree = 4
tolerence = 20
# directories of SPM results and VAF files
addrres = "../DataPreprocessing/data/SPM_out/"
addrvaf = "../DataPreprocessing/data/COMB/"
#--------------------------

res = read.table(paste0(addrres, "SPM_ccs_out_K", nptwy, "_T", ntree, "_Tol", tolerence, "_", in_filename), sep = '\t')
res = as.matrix(res)
vaf = read.table(paste0(addrvaf, in_filename), sep = '\t')
m = dim(vaf)[1]
n = dim(vaf)[2]

base_s = n*nptwy +m*n*nptwy
base_a = n*nptwy
base_A = n*nptwy +m*n*nptwy + m*ntree
base_p = 0
  
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
  }
  
  for (t in 1:ntree){
    part = c()
    for (i in 1:m){
      if (res[r, base_s+ (i-1)*ntree + t] == 1)
        part = c(part, i)
    }
    print(paste(c("subtype", t, "consists of", length(part), "samples:", part), collapse = " "))
    print(paste(c("subtype", t, "is in the shape"), collapse = " "))
    for (k in 1:nptwy) {
      str = base_A+ (t-1)*(nptwy+1)*(nptwy+1) + (k-1)*(nptwy+1) + 1
      end = base_A+ (t-1)*(nptwy+1)*(nptwy+1) + k*(nptwy+1) -1
      print(unname(res[r,str:end]))
    }
  }
  print("-------------------")
}
