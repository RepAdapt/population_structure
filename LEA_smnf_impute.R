####################
#
#
#
# Outline
# -------
# 
# 
# TODO
# ----
# get imputed back to VCF
#
####################

library(LEA)
library(bigsnpr)

print(sessionInfo())

# get input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: LEA_smnf_impute.R <input.vcf> <impute_dir>")
}
vcf <- args[1]
impute_dir <- args[2]

cat(sprintf('vcf = %s', vcf))

# set global vars
setwd(impute_dir)
name <- sub("\\.(vcf|vcf\\.gz)$", "", basename(vcf), perl = TRUE)

# convert vcf to geno file, then to lfmm
cat(sprintf('\nconverting vcf to lfmm %s', Sys.time()))
geno_file = sprintf("%s.geno", name)
lfmm_file = sprintf("%s.lfmm", name)

vcf2geno(vcf, geno_file)
geno2lfmm(geno_file, lfmm_file)


# run snmf
choose_k <- function(LFMM_file, Kvals=3:3, reps=2, project="new"){  # UPDATE KVALS, REPS !!!!!!!!!!
    ################################################################################################
    # 
    # Estimate individual ancestry coefficients and ancestral allele frequencies.
    #
    # Parameters
    # ----------
    # lfmm_file : path
    #    path to data in .lfmm format
    # Kvals : vector
    #    vector of k values to evaluate
    # reps : int
    #    An integer corresponding with the number of repetitions for each value of ‘K’.
    # project : str
    #    create (or overwrite) a new snmf project
    # 
    # Returns
    # -------
    # list 
    #    a list with the following 1) snmfProject class object, bestK (int), best_run_idx (int)
    ################################################################################################
    cat(sprintf('\nrunning snmf %s\n', Sys.time()))
    proj <- snmf(LFMM_file, K = Kvals, repetitions = reps, entropy = TRUE, project = project)
    
    # choose K
    cat(sprintf('\nchoosing K %s\n', Sys.time()))
    ce_mat <- sapply(Kvals, function(k) cross.entropy(proj, K = k))
    ce_med <- apply(ce_mat, 2, median)
    
    bestK <- Kvals[which.min(ce_med)]
    run_ce <- cross.entropy(proj, K = bestK)
    best_run_idx <- which.min(run_ce)
    
    cat(sprintf("best K = %s (run %s)\n", bestK, best_run_idx))

    return(list(proj = proj, bestK = bestK, best_run_idx = best_run_idx))
}

proj_and_best_k_results <- choose_k(lfmm_file, project='new')


# impute
cat(sprintf('\nimputing %s\n', Sys.time()))
x <- impute(  # returns NULL
    object = proj_and_best_k_results$proj,
    K = proj_and_best_k_results$bestK,
    run = proj_and_best_k_results$best_run_idx,
    input.file = lfmm_file
)
# rename output
imputed_file = sprintf("%s_imputed.lfmm", name)
cat(sprintf('\nrenaming imputed output file: %s\n', imputed_file))
file.rename(
    sprintf("%s.lfmm_imputed.lfmm", name),
    imputed_file
)

# BIGSNPR
cat(sprintf('\nstarting bigsnpr\n'))
X <- read.lfmm(imputed_file)
stopifnot(sum(is.na(X)) == 0)
stopifnot(all(X %in% 0:2))

bk <- sprintf("%s_imputed", name)
if (file.exists(bk)){
    # avoid error about an existing backing file
    file.remove(bk)
}
G <- bigstatsr::FBM.code256(nrow(X), ncol(X), backingfile = bk, code = CODE_IMPUTE_PRED)
G[] <- X

fam <- read.table(sprintf('%s.fam', name))
bim <- read.table(sprintf("%s.bim", name), stringsAsFactors = FALSE)
                 
big_svd <- snp_autoSVD(
    G,
    infos.chr=bim[[1]],
    infos.pos=bim[[4]],
    thr.r2=0.2,    # UPDATE!
    size=100/0.2,  # UPDATE!!
    min.maf=0.01,  # MAY BE REDUNDANT IF I FILTER FOR MAF ABOVE
    roll.size=10   # UPDATE!!!!!!
)

# compute PCs
pcs <- predict(big_svd)
colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
row.names(pcs) <- fam[, 'V1']  # add sample names to rows of PC loadings
pcs_file <- sprintf("%s_imputed_pcs.txt", name)
write.table(pcs, pcs_file, sep = "\t")
cat(sprintf('\nwrote PC loadings to : %s\n', pcs_file))


# ADMIXTURE ANALYSIS




# remove unnecessary files
file.remove(geno_file)
file.remove(lfmm_file)