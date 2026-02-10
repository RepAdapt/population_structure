####################
#
# Impute and thin filtered SNPs from an input VCF file.
#
# Usage
# -----
# Rscript LEA_smnf_impute.R vcf_path output_path cpus
#
# Parameters
# ----------
# vcf : Path
#     path to unzipped filtered vcf (basename.vcf) - the basename will be used to prefix output files
#     vcf is assumed to be filtered according to PopStruct_01_filter_vcf.py
#         - remove individuals with > 10% missing data across SNPs, then
#         - remove SNPs with > 10% missing data across individuals
# output_path : Path
#     Where to save output
# threads : int
#     the number of cpus for parallelization of snmf function
# 
# Notes
# -----
# - if run manually, direct output to a file : Rscript LEA_smnf_impute.R vcf_path output_path cpus > output
# - after running, cat either the slurm.out file or manually created output file for useful information:
#     cat output | grep INFO
# - admixture figures are available for before and after thinning
# 
# Assumed
# -------
# .bim, .bam, .fam files created in output_path from:
#    `plink --vcf basename.vcf --make-bed --out basename --keep-allele-order --allow-extra-chr`
#    in this pipeline, the plink command this is run upstream (see PopStruct_02_impute_filter_thin.py)
#
# Outline
# -------
# 1. GET INPUT ARGUMENTS, METADATA, AND FUNCTIONS
# 2. CONVERT VCF TO GENO FILE, THEN TO LFMM
# 3. RUN SNMF - PLOT ANCESTRY PROPORTIONS
# 4. IMPUTE, FILTER MAF (FULL SNP DATASET)
# 5. BIGSNPR
# 6. ESTIMATE PCs FOR THINNED SNP DATA
# 7. ISOLATE THINNED SNP SET - PLOT AND SAVE ANCESTRY PROPORTIONS
# 8. ADMIXTURE ANALYSIS ON THINNED SNP DATASET
#
####################

library(LEA)
library(bigsnpr)

print(sessionInfo())

# 1. GET INPUT ARGUMENTS, METADATA, AND FUNCTIONS
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage: LEA_smnf_impute.R <input.vcf> <output_path> <cpus>")
}
vcf <- args[1]
output_path <- args[2]
cpus <- as.integer(args[3])
setwd(output_path)

source(sprintf('%s/pop_struct/LEA_config.R', Sys.getenv('HOME'))  # args for snp_autoSVD

cat(sprintf('\nINFO [%s] vcf = %s\n', Sys.time(), vcf))

## set global vars
name <- sub("\\.(vcf|vcf\\.gz)$", "", basename(vcf), perl = TRUE)

## metadata
fam <- read.table(sprintf('%s.fam', name))
bim <- read.table(sprintf("%s.bim", name), stringsAsFactors = FALSE)

## functions
choose_k <- function(LFMM_file, Kvals=1:10, reps=10, project="new", CPU=1){
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
    #    if 'new' - current project is removed and a new one is created to store the result
    # 
    # Returns
    # -------
    # list 
    #    a list with the following 1) snmfProject class object, bestK (int), best_run_idx (int)
    ################################################################################################
    proj <- snmf(LFMM_file, K = Kvals, repetitions = reps, entropy = TRUE, project = project, CPU=CPU)
    
    # choose K
    cat(sprintf('\nINFO [%s] choosing K\n', Sys.time()))
    ce_mat <- sapply(Kvals, function(k) cross.entropy(proj, K = k))
    ce_med <- apply(ce_mat, 2, median)
    
    bestK <- Kvals[which.min(ce_med)]
    run_ce <- cross.entropy(proj, K = bestK)
    best_run_idx <- which.min(run_ce)
    
    cat(sprintf("INFO [%s] best K = %s (run %s)\n", Sys.time(), bestK, best_run_idx))

    return(list(proj = proj, bestK = bestK, best_run_idx = best_run_idx))
}

make_ancestry_plot <- function(file_basename, proj, K, run){
    #################################################################
    #
    # Create a color blind-friendly ancestry plot from snmf output
    #
    # Parameters
    # ----------
    # filename : /path/to/save.pdf
    # proj : snmfProject class object
    # K : best K determined from cross-entropy
    # run : best run index where the best K was determined
    #################################################################
    colors <- c("#56B4E9", "#ECE237", "#949494", "#FBAFF6", "#CA9161", "#CC78BC",
                "#D55E00", "#029E73", "#DE8F05", "#0373B2", "black", "white")

    filename <- sprintf("%s.pdf", file_basename)
    pdf(filename, height=4, width=8)
    barchart(
        proj,
        K = K,
        run = run,
        border = NA,
        space = 0,
        col = colors,
        xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = sprintf("Ancestry matrix\n%s", file_basename)
    ) -> bp
    dev.off()
    cat(sprintf('\nINFO [%s] ancestry pdf saved to: %s\n', Sys.time(), filename))
}

                     
# 2. CONVERT VCF TO GENO FILE, THEN TO LFMM
cat(sprintf('\nINFO [%s] converting vcf to lfmm\n', Sys.time()))
geno_file = sprintf("%s.geno", name)
lfmm_file = sprintf("%s.lfmm", name)

vcf2geno(vcf, geno_file)
geno2lfmm(geno_file, lfmm_file)


# 3. RUN SNMF - PLOT ANCESTRY PROPORTIONS
cat(sprintf('\nINFO [%s] running snmf on unimputed snps\n', Sys.time()))
proj_and_best_k_results <- choose_k(lfmm_file, project='new', CPU=cpus)
make_ancestry_plot(
    file_basename = sprintf('%s_ancestry_unimputed', name),
    proj = proj_and_best_k_results$proj,
    K = proj_and_best_k_results$bestK,
    run = proj_and_best_k_results$best_run_idx
)


# 4. IMPUTE, FILTER MAF (FULL SNP DATASET)
cat(sprintf('\nINFO [%s] imputing snps\n', Sys.time()))
x <- impute(  # returns NULL
    object = proj_and_best_k_results$proj,
    K = proj_and_best_k_results$bestK,
    run = proj_and_best_k_results$best_run_idx,
    input.file = lfmm_file
)

## rename output
imputed_file = sprintf("%s_imputed.lfmm", name)
cat(sprintf('\nINFO [%s] renaming imputed output file: %s\n', Sys.time(), imputed_file))
file.rename(
    sprintf("%s.lfmm_imputed.lfmm", name),
    imputed_file
)

## filter for MAF and save
X <- read.lfmm(imputed_file)
row.names(X) <- fam[, 'V1']    # add sample names to rows of X
colnames(X) <- bim[, 'V2']     # add locus names to cols of X
stopifnot(sum(is.na(X)) == 0)  # check nothing is NA
stopifnot(all(X %in% 0:2))     # check good genotypes

af <- colMeans(X, na.rm = TRUE) / 2
keep_loci <- af >= min_maf & af <= (1 - min_maf)  # min_maf sourced from LEA_config.R
X_filtered <- X[, keep_loci]
rm(X)  # avoid accidental re-use
filtered_file <- sprintf("%s_imputed_maf-filtered.txt", name)
write.table(X_filtered, filtered_file, sep='\t', row.names=TRUE, col.names=TRUE)
cat(sprintf('\nINFO [%s] wrote full snp data file to: %s\n', Sys.time(), filtered_file))


# 5. BIGSNPR
cat(sprintf('\nINFO [%s] starting bigsnpr\n', Sys.time()))

## create big snp class object
bk <- sprintf("%s_imputed_maf-filtered", name)
if (file.exists(sprintf('%s.bk', bk))){
    # avoid error about an existing backing file
    null <- file.remove(sprintf('%s.bk', bk))
}
G <- bigstatsr::FBM.code256(nrow(X_filtered), ncol(X_filtered), backingfile = bk, code = CODE_IMPUTE_PRED)
G[] <- X_filtered

## run snp_autoSVD on imputed data
cat(sprintf('\nINFO [%s] running snp_autoSVD on full snp data\n', Sys.time()))
filtered_bim = bim[bim[[2]] %in% names(keep_loci)[keep_loci], ]
big_svd <- snp_autoSVD(
    G,
    infos.chr=filtered_bim[[1]],
    infos.pos=filtered_bim[[4]],
    thr.r2=thr_r2,        # sourced from LEA_config.R
    size=size,            # sourced from LEA_config.R
    min.maf=min_maf,      # sourced from LEA_config.R - should be redundant because of upstream filtering
    # roll.size=roll_size,  # sourced from LEA_config.R - for debugging
)


# 6. ESTIMATE PCs FOR THINNED SNP DATA
pcs <- predict(big_svd)
colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
row.names(pcs) <- fam[, 'V1']  # add sample names to rows of PC loadings
                     
pcs_file <- sprintf("%s_imputed_maf-filtered_thinned_pcs.txt", name)
write.table(pcs, pcs_file, sep = "\t")
cat(sprintf('\nINFO [%s] wrote PC loadings of thinned snp data to : %s\n', Sys.time(), pcs_file))


# 7. ISOLATE THINNED SNP SET - PLOT AND SAVE ANCESTRY PROPORTIONS
thinned_indices = attr(big_svd, 'subset')
thinned_snps = X_filtered[, thinned_indices]
                     
thinned_file = sprintf("%s_imputed_maf-filtered_thinned.txt", name)
write.table(thinned_snps, thinned_file, sep='\t', row.names=TRUE, col.names=TRUE)
                    
thinned_lfmm = sprintf("%s_imputed_maf-filtered_thinned.lfmm", name)
write.lfmm(thinned_snps, thinned_lfmm)
cat(sprintf('\nINFO [%s] wrote thinned snp data to : %s\n', Sys.time(), thinned_file))


# 8. ADMIXTURE ANALYSIS ON THIN SNP DATASET - PLOT AND SAVE ANCESTRY PROPORTIONS
cat(sprintf('\nINFO [%s] running smnf on thinned imputed snps\n', Sys.time()))

thinimp_proj_and_best_k_results <- choose_k(lfmm_file, Kvals=proj_and_best_k_results$bestK, project='new', CPU=cpus)

## plot ancestry proportions
# ancestry_pdf <- sprintf("%s_imputed_maf-filtered_thinned_ancestry.pdf", name)
                 
make_ancestry_plot(
    file_basename = sprintf("%s_imputed_maf-filtered_thinned_ancestry", name),
    proj = thinimp_proj_and_best_k_results$proj,
    K = proj_and_best_k_results$bestK,  # use K from unimputed SNPs
    run = thinimp_proj_and_best_k_results$best_run_idx
)

## save final ancestry proportions
anc_props <- Q(
    thinimp_proj_and_best_k_results$proj,
    K = proj_and_best_k_results$bestK,  # use K from unimputed SNPs
    run = thinimp_proj_and_best_k_results$best_run_idx
)
row.names(anc_props) <- fam[, 'V1']  # add sample names to rows of PC loadings

props_file <- sprintf("%s_imputed_maf-filtered_thinned_ancestry_proportions.txt", name)
write.table(anc_props, props_file, sep='\t', row.names=T, col.names=T)
cat(sprintf('\nINFO [%s] ancestry proportions saved to: %s\n', Sys.time(), props_file))

                     
print('DONE!')
