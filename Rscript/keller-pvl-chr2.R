# The code below performs the two-dimensional QTL scan over the specified genomic region. 
# We analyze the chromosome 2 hotspot from Keller et al 2018 (GENETICS).
# Each scan considers two traits. One of those two is the Hnf4a gene expression level. 
# Keller et al 2018 state that Hnf4a is causal for 88 of the 147 hotspot nonlocal traits.

##First read in the arguments listed at the command line

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
print(args)
##args is now a list of character vectors
print(args$argname)
proc_num <- as.numeric(args$argname)
print(proc_num)
(hot_indic <- proc_num %% 139 + 1) # define hot_indic 
(local_indic <- proc_num %/% 139 + 1) # define local_indic
run_num <- as.numeric(args$run_num)
print(run_num)
(nsnp <- as.numeric(args$nsnp))
(s1 <- as.numeric(args$s1))

###############


# load expression traits
readRDS("data/hotspot_expr_tib2_keller_chr2.rds") -> hotspot
readRDS("data/keller2018-chr2-local-expr.rds") -> local

# load chr2 allele probabilities
readRDS("genoprobs_chr2.rds") -> geno # genoprobs_chr2.rds is on SQUID

# load kinship matrix (LOCO, ie, for chromosome 2, ie, doesn't use chr2 data)
readRDS("data/kinship_chr2.rds") -> kinship

# load covariates
readRDS("data/addcovar.rds") -> covar


# remove subjects with missing data

id2keep <- local$mouse_id
gg <- geno[[1]]
gg2 <- gg[rownames(gg) %in% id2keep, , ]
kk <- kinship[[1]]
kk2 <- kk[rownames(kk) %in% id2keep, colnames(kk) %in% id2keep]
cc2 <- covar[rownames(covar) %in% id2keep, ]


# create matrix of two expression traits
hnf4a_id <- "ENSMUSG00000017950"
local2 <- local[ , !colnames(local) == hnf4a_id]
pheno <- cbind(local2[ , local_indic, drop = FALSE], hotspot[ , hot_indic, drop = FALSE])
rownames(pheno) <- local$mouse_id
# verify that names match in all objects
sum(rownames(pheno) == rownames(gg2))
sum(rownames(pheno) == rownames(kk2))
sum(rownames(pheno) == colnames(kk2))
sum(rownames(pheno) == rownames(cc2))
phenames <- c(colnames(local2)[local_indic], colnames(hotspot)[hot_indic])
# two-dimensional scan

library(qtl2pleio)
s_out <- scan_pvl(probs = gg2,
         pheno = pheno,
         kinship = kk2,
         covariates = cc2[ , -5], # need to remove column 5 because we have no mice from wave 5
         start_snp1 = s1,
         n_snp = nsnp
         #n_snp = 5
           )

colnames(pheno)
fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
write.table(s_out, fn_out, quote = FALSE)
q("no")
