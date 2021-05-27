library(data.table)
library(foreach)
library(doMC)
library(ggplot2)
library(ggrastr)
library(readstata13)

options(ggrastr.default.dpi=1200)

# Options for plink:
# WARNING: this takes a long time to run > 12 hours due to the size of the UKB genotype data:
# You either need to extract the variants of interest prior to GWAS (slow I/O) or 
# plink will calculate allele frequencies for all variants (also slow).
nCores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
memory <- as.integer(nCores) * 6000 # MB

registerDoMC(nCores)
setDTthreads(nCores)

# Create output directories
if (!dir.exists("paper_output")) dir.create("paper_output")
if (!dir.exists("gwas")) dir.create("gwas")

# Load data with adjusted values at each step
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Load technical information and filter to samples in UKB raw data
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(sample_id, visit)], on = .(sample_id, visit)]

# For samples with data at both baseline and repeat assessment, take baseline data
baseline <- sinfo[!(sample_removed) & (in_ukb_raw) & visit == "Main Phase"]
repeat_visit <- sinfo[!(sample_removed) & (in_ukb_raw) & visit == "Repeat Assessment"]
repeat_visit <- repeat_visit[!baseline, on = .(eid_7439)]
first <- rbind(baseline, repeat_visit)
dat <- dat[first[,.(eid_7439, sample_id, visit)], on = .(eid_7439, sample_id, visit)]

# Load all the data needed for GWAS
geno_pcs <- fread("data/genetic_reference/ukb_sqc_v2.txt")
geno_pcs <- geno_pcs[,.(eid, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

kinship <- fread("data/genetic_reference/kinship_relatedness.txt")

age <- read.dta13("data/ceu_curated_phenotypes/20210302/STATA/repeats.dta")
age <- as.data.table(age)
age <- age[,.(eid_7439=idno, age=ages_rep, visit=fcase(repno == 0, "Main Phase", repno == 1, "Repeat Assessment", default = NA))]
age <- age[!is.na(visit)]

pheno <- read.dta13("data/ceu_curated_phenotypes/20210302/STATA/analysis.dta")
pheno <- as.data.table(pheno)
pheno <- pheno[,.(eid_7439=idno, white_british=racebin == "White")]

psam <- fread("data/imputed_pgen/ukb_imp_v3_dedup_chr22.psam")

# Build table of phenotype data
pheno <- pheno[age, on = .(eid_7439), nomatch=0]
pheno[,eid_7439 := as.integer(eid_7439)]
pheno <- pheno[geno_pcs, on = .(eid_7439=eid), nomatch=0]
pheno <- pheno[first[,.(eid_7439, visit)], on = .(eid_7439, visit), nomatch=0]

pheno[dat[variable == "Ala"], on = .(eid_7439, visit), raw_ala := log(raw)]
pheno[dat[variable == "His"], on = .(eid_7439, visit), raw_his := log(raw)]
pheno[dat[variable == "Gln"], on = .(eid_7439, visit), raw_gln := log(raw)]
pheno[dat[variable == "Albumin"], on = .(eid_7439, visit), raw_alb := log(raw)]

pheno[dat[variable == "Ala"], on = .(eid_7439, visit), postqc_ala := log(adj5)]
pheno[dat[variable == "His"], on = .(eid_7439, visit), postqc_his := log(adj5)]
pheno[dat[variable == "Gln"], on = .(eid_7439, visit), postqc_gln := log(adj5)]
pheno[dat[variable == "Albumin"], on = .(eid_7439, visit), postqc_alb := log(adj5)]

# Filter to participants of White British ancestry
pheno <- pheno[(white_british)]

# From pairs of close relatives (first- or second- degree), pick one from each pair
kinship <- kinship[ID1 %in% pheno$eid_7439 & ID2 %in% pheno$eid_7439]
pheno <- pheno[!kinship[Kinship > 0.0884], on = .(eid_7439=ID2)] # cutoff from KING manual http://people.virginia.edu/~wc9c/KING/manual.html

# Write out global sample keep file for genetic data extraction
keep <- pheno[, .(`#FID`=eid_7439, IID=eid_7439)]
fwrite(keep, sep="\t", quote=FALSE, file="gwas/samples_passing_qc.txt")

# Write out phenotype data and sample keep files
ala_raw_keep <- pheno[!is.na(raw_ala), .(`#FID`=eid_7439, IID=eid_7439)]
ala_raw_pheno <- pheno[!is.na(raw_ala), .(`#FID`=eid_7439, IID=eid_7439, raw_ala)]
ala_raw_covar <- pheno[!is.na(raw_ala), .(`#FID`=eid_7439, IID=eid_7439, age, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

fwrite(ala_raw_keep, sep="\t", quote=FALSE, file="gwas/ala_raw_keep.txt")
fwrite(ala_raw_pheno, sep="\t", quote=FALSE, file="gwas/ala_raw_pheno.txt")
fwrite(ala_raw_covar, sep="\t", quote=FALSE, file="gwas/ala_raw_covar.txt")

his_raw_keep <- pheno[!is.na(raw_his), .(`#FID`=eid_7439, IID=eid_7439)]
his_raw_pheno <- pheno[!is.na(raw_his), .(`#FID`=eid_7439, IID=eid_7439, raw_his)]
his_raw_covar <- pheno[!is.na(raw_his), .(`#FID`=eid_7439, IID=eid_7439, age, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

fwrite(his_raw_keep, sep="\t", quote=FALSE, file="gwas/his_raw_keep.txt")
fwrite(his_raw_pheno, sep="\t", quote=FALSE, file="gwas/his_raw_pheno.txt")
fwrite(his_raw_covar, sep="\t", quote=FALSE, file="gwas/his_raw_covar.txt")

gln_raw_keep <- pheno[!is.na(raw_gln), .(`#FID`=eid_7439, IID=eid_7439)]
gln_raw_pheno <- pheno[!is.na(raw_gln), .(`#FID`=eid_7439, IID=eid_7439, raw_gln)]
gln_raw_covar <- pheno[!is.na(raw_gln), .(`#FID`=eid_7439, IID=eid_7439, age, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

fwrite(gln_raw_keep, sep="\t", quote=FALSE, file="gwas/gln_raw_keep.txt")
fwrite(gln_raw_pheno, sep="\t", quote=FALSE, file="gwas/gln_raw_pheno.txt")
fwrite(gln_raw_covar, sep="\t", quote=FALSE, file="gwas/gln_raw_covar.txt")

alb_raw_keep <- pheno[!is.na(raw_alb), .(`#FID`=eid_7439, IID=eid_7439)]
alb_raw_pheno <- pheno[!is.na(raw_alb), .(`#FID`=eid_7439, IID=eid_7439, raw_alb)]
alb_raw_covar <- pheno[!is.na(raw_alb), .(`#FID`=eid_7439, IID=eid_7439, age, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

fwrite(alb_raw_keep, sep="\t", quote=FALSE, file="gwas/alb_raw_keep.txt")
fwrite(alb_raw_pheno, sep="\t", quote=FALSE, file="gwas/alb_raw_pheno.txt")
fwrite(alb_raw_covar, sep="\t", quote=FALSE, file="gwas/alb_raw_covar.txt")

ala_postqc_keep <- pheno[!is.na(postqc_ala), .(`#FID`=eid_7439, IID=eid_7439)]
ala_postqc_pheno <- pheno[!is.na(postqc_ala), .(`#FID`=eid_7439, IID=eid_7439, postqc_ala)]
ala_postqc_covar <- pheno[!is.na(postqc_ala), .(`#FID`=eid_7439, IID=eid_7439, age, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

fwrite(ala_postqc_keep, sep="\t", quote=FALSE, file="gwas/ala_postqc_keep.txt")
fwrite(ala_postqc_pheno, sep="\t", quote=FALSE, file="gwas/ala_postqc_pheno.txt")
fwrite(ala_postqc_covar, sep="\t", quote=FALSE, file="gwas/ala_postqc_covar.txt")

his_postqc_keep <- pheno[!is.na(postqc_his), .(`#FID`=eid_7439, IID=eid_7439)]
his_postqc_pheno <- pheno[!is.na(postqc_his), .(`#FID`=eid_7439, IID=eid_7439, postqc_his)]
his_postqc_covar <- pheno[!is.na(postqc_his), .(`#FID`=eid_7439, IID=eid_7439, age, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

fwrite(his_postqc_keep, sep="\t", quote=FALSE, file="gwas/his_postqc_keep.txt")
fwrite(his_postqc_pheno, sep="\t", quote=FALSE, file="gwas/his_postqc_pheno.txt")
fwrite(his_postqc_covar, sep="\t", quote=FALSE, file="gwas/his_postqc_covar.txt")

gln_postqc_keep <- pheno[!is.na(postqc_gln), .(`#FID`=eid_7439, IID=eid_7439)]
gln_postqc_pheno <- pheno[!is.na(postqc_gln), .(`#FID`=eid_7439, IID=eid_7439, postqc_gln)]
gln_postqc_covar <- pheno[!is.na(postqc_gln), .(`#FID`=eid_7439, IID=eid_7439, age, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

fwrite(gln_postqc_keep, sep="\t", quote=FALSE, file="gwas/gln_postqc_keep.txt")
fwrite(gln_postqc_pheno, sep="\t", quote=FALSE, file="gwas/gln_postqc_pheno.txt")
fwrite(gln_postqc_covar, sep="\t", quote=FALSE, file="gwas/gln_postqc_covar.txt")

alb_postqc_keep <- pheno[!is.na(postqc_alb), .(`#FID`=eid_7439, IID=eid_7439)]
alb_postqc_pheno <- pheno[!is.na(postqc_alb), .(`#FID`=eid_7439, IID=eid_7439, postqc_alb)]
alb_postqc_covar <- pheno[!is.na(postqc_alb), .(`#FID`=eid_7439, IID=eid_7439, age, genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

fwrite(alb_postqc_keep, sep="\t", quote=FALSE, file="gwas/alb_postqc_keep.txt")
fwrite(alb_postqc_pheno, sep="\t", quote=FALSE, file="gwas/alb_postqc_pheno.txt")
fwrite(alb_postqc_covar, sep="\t", quote=FALSE, file="gwas/alb_postqc_covar.txt")

# Filter to bi-allelic SNPs with MAF > 1% and info > 0.4
for (chr in 1:22) {
  info <- fread(cmd=sprintf("grep -v '^#' data/genetic_reference/ukb_impv3_chr%s_snpstats.txt", chr))
  # Filter on individual allele frequencies instead of MAF so we don't through out bi-allelic SNPs with rare variants that make them multi-allelic
  info <- info[nchar(alleleA) == 1 & nchar(alleleB) == 1 & alleleA_frequency > 0.01 & alleleB_frequency > 0.01 & impute_info > 0.4]
  multi <- info[,.N,by=rsid][N > 1]
  info <- info[!multi, on = .(rsid)]

  pvar <- fread(sprintf("data/imputed_pgen/ukb_imp_v3_dedup_chr%s.pvar", chr))
  info <- rbind(
    info[pvar, on = .(chromosome=`#CHROM`, position=POS, alleleA=ALT, alleleB=REF), nomatch=0],
    info[pvar, on = .(chromosome=`#CHROM`, position=POS, alleleA=REF, alleleB=ALT), nomatch=0]
  )
  
  fwrite(info[,.(ID)], quote=FALSE, col.names=FALSE, file=sprintf("gwas/chr%s_keep.txt", chr))
}

# Extract variant and sample subsets to dramatically speed up the GWAS by skipping allele frequency
# calculations for excluded variants.
foreach(chr = 22:1) %dopar% {
  cmd <- sprintf("plink2 --pfile data/imputed_pgen/ukb_imp_v3_dedup_chr%s", chr)
  cmd <- sprintf("%s --extract gwas/chr%s_keep.txt", cmd, chr)
  cmd <- sprintf("%s --keep gwas/samples_passing_qc.txt", cmd)
  cmd <- sprintf("%s --make-pgen --out gwas/pass_qc_chr%s", cmd, chr)
	cmd <- sprintf("%s --threads %s --memory %s", cmd, 1, memory / nCores) # memory and cores set at top of file

	# Run and wait
	system(cmd, wait=TRUE)
}

# Run association scan in plink
for (pheno in c("ala", "his", "gln", "alb")) {
  for (type in c("raw", "postqc")) {
		for (chr in 22:1) {
			cmd <- sprintf("plink2 --pfile gwas/pass_qc_chr%s", chr)
			cmd <- sprintf("%s --maf 0.01", cmd) # re-apply MAF filter in the subset of samples with biomarker data passing QC.
			cmd <- sprintf("%s --pheno gwas/%s_%s_pheno.txt", cmd, pheno, type)
			cmd <- sprintf("%s --covar gwas/%s_%s_covar.txt", cmd, pheno, type)
      cmd <- sprintf("%s --keep gwas/%s_%s_keep.txt", cmd, pheno, type)
			cmd <- sprintf("%s --glm sex log10 hide-covar", cmd)
			cmd <- sprintf("%s cols=chrom,pos,ref,alt,a1freq,test,nobs,beta,se,p,err", cmd)
			cmd <- sprintf("%s --quantile-normalize --pfilter 0.001", cmd)
			cmd <- sprintf("%s --out gwas/%s_%s_chr%s", cmd, pheno, type, chr)
			cmd <- sprintf("%s --threads %s --memory %s", cmd, nCores, memory) # System resource parameters set at top of file
		 
			# Run and wait
			system(cmd, wait=TRUE)
		}
  }
}

# Collate plink logs
system("touch paper_output/gwas_plink_logs.txt", wait=TRUE)
for (chr in 1:22) {
  system(sprintf("cat gwas/pass_qc_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)

  system(sprintf("cat gwas/ala_raw_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)

  system(sprintf("cat gwas/ala_postqc_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)

  system(sprintf("cat gwas/his_raw_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)

  system(sprintf("cat gwas/his_postqc_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)

  system(sprintf("cat gwas/gln_raw_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)

  system(sprintf("cat gwas/gln_postqc_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)

  system(sprintf("cat gwas/alb_raw_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)

  system(sprintf("cat gwas/alb_postqc_chr%s.log >> paper_output/gwas_plink_logs.txt", chr), wait=TRUE)
  system("echo '' >> paper_output/gwas_plink_logs.txt", wait=TRUE)
}

# Extract results
gwas <- foreach(pheno = c("ala", "his", "gln", "alb"), .combine=rbind) %do% {
  foreach(type = c("raw", "postqc"), .combine=rbind) %do% {
    collated <- foreach(chr = 1:22, .combine=rbind) %do% {
      fread(sprintf("gwas/%s_%s_chr%s.%s_%s.glm.linear", pheno, type, chr, type, pheno))
    }
    fwrite(collated, sep="\t", quote=FALSE, file=sprintf("paper_output/gwas_results_%s_%s.txt", type, pheno))
    cbind(variable = pheno, type = type, collated)
  }
}

# No longer need gwas directory
system("rm -rf gwas", wait=TRUE)

# Compare P-values from raw vs. postqc
pcomp <- dcast(gwas, variable + ID ~ type, value.var="LOG10_P")
pcomp <- pcomp[!is.na(raw) & !is.na(postqc)] # must be P < 0.001 for both raw and postqc concentrations to compare

g <- ggplot(pcomp) +
  aes(x=raw, y=postqc) +
  rasterize(geom_point(shape=21, color="white", fill="#304830", stroke=0.1, size=1)) +
  geom_abline(intercept=0, slope=1, linetype=2, color="black", size=0.3) +
  xlab("Raw -log10 P-values") +
  ylab("Post-QC -log10 P-values") +
  facet_wrap(~ variable, scales="free", nrow=1) +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(), strip.text=element_text(size=7)
  )

ggsave(g, width=7.2, height=2, file="paper_output/gwas_compare.pdf")

# For Manhattan plots, we need to set a cumulative genome location, for which we need
# to know the size of each chromosome
chr_info <- foreach(chr = 1:22, .combine=rbind) %do% {
  prefix <- "data/imputed_pgen/ukb_imp_v3_dedup_chr"
  last <- fread(cmd = sprintf("tail -1 %s%s.pvar | cut -f 1,2", prefix, chr))
  setnames(last, c("chr", "last_pos"))
  last
}

chr_info[, cumul_chr_start := c(0, cumsum(as.numeric(last_pos))[-22])]

gwas[chr_info, on = .(`#CHROM` = chr), cumul_pos := POS + cumul_chr_start]

# Make manhattan plots
g <- ggplot(gwas[type == "raw"]) + 
  aes(x=cumul_pos, y=LOG10_P, color=factor(`#CHROM` %% 2)) +
  rasterize(geom_point(shape=19, stroke=0, size=1)) +
  geom_hline(yintercept=-log10(5e-8), color="red", linetype=2, size=0.4) +
  scale_x_continuous(expand=c(0,0), limits=c(1, chr_info[22, cumul_chr_start + last_pos])) +
  scale_y_continuous(name="-log10 P-value", limits=c(0, NA), expand = expansion(mult = c(0, .1))) +
  facet_wrap(~ variable, ncol=1, scales="free") +
  scale_color_manual(values=c("0"="#969696", "1"="#fd8d3c")) +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(), strip.text=element_text(size=7),
        legend.position="none"
  )
ggsave(g, width=7.2, height=7.2, file="paper_output/raw_gwas_manhattans.pdf")

g <- ggplot(gwas[type == "postqc"]) + 
  aes(x=cumul_pos, y=LOG10_P, color=factor(`#CHROM` %% 2)) +
  rasterize(geom_point(shape=19, stroke=0, size=1)) +
  geom_hline(yintercept=-log10(5e-8), color="red", linetype=2, size=0.4) +
  scale_x_continuous(expand=c(0,0), limits=c(1, chr_info[22, cumul_chr_start + last_pos])) +
  scale_y_continuous(name="-log10 P-value", limits=c(0, NA), expand = expansion(mult = c(0, .1))) +
  facet_wrap(~ variable, ncol=1, scales="free") +
  scale_color_manual(values=c("0"="#969696", "1"="#fd8d3c")) +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(), strip.text=element_text(size=7),
        legend.position="none"
  )
ggsave(g, width=7.2, height=7.2, file="paper_output/postqc_gwas_manhattans.pdf")

# Load in approximately independent LD blocks to help with peak annotation
ld_blocks <- fread("data/Berisa_2016_LD_blocks/EUR_1000G_ind_ld_blocks.bed")
ld_blocks[, chr := as.integer(gsub("chr", "", chr))]
ld_blocks[, block_num := .I]
gwas[ld_blocks, on = .(`#CHROM`=chr, POS >= start, POS < stop), ld_block := i.block_num]

peaks <- gwas[LOG10_P > -log10(5e-8), .SD[which.max(LOG10_P)], by=.(variable, type, ld_block)]

# Query Variant Effect Predictor for each peak
library(httr)
library(jsonlite)
library(xml2)

r <- POST("https://grch37.rest.ensembl.org/vep/human/id", content_type("application/json"), accept("application/json"), 
          body = sprintf('{ "ids" : [ %s ] }', paste(sprintf('"%s"', unique(peaks$ID)), collapse=",")))

stop_for_status(r)

# Extract
unlist_json <- function(l) {
	foreach(li = seq_along(l), .combine=c) %do% {
    if (is.null(l[[li]])) {
      return(NA)
    } else if (!is.null(dim(l[[li]]))) {
      return(paste(as.vector(l[[li]]), collapse=";"))
    } else {
      return(l[[li]])
    }
  }
}

rbindf <- function(...) { rbind(..., fill=TRUE) }
vep <- foreach(ii = seq_along(content(r)), .combine=rbindf) %do% {
  rsID <- content(r)[[ii]][["id"]]
  most_severe <- content(r)[[ii]][["most_severe_consequence"]]
  consequences <- content(r)[[ii]][["transcript_consequences"]]
  consequences <- fromJSON(toJSON(consequences))
  consequences <- as.data.table(lapply(consequences, unlist_json))
  cbind(ID = rsID, most_severe = most_severe, consequences)
}

# Filter to most severe consequences per ID
vep <- vep[, .SD[consequence_terms %like% unique(most_severe)], by=ID]

# Extract columns and rows of interest
vep <- vep[, .(ID, Gene=gene_symbol, biotype, distance, most_severe, consequence_terms, variant_allele, impact, polyphen_prediction, sift_prediction)]
vep <- unique(vep) # drop duplicates arising from different ensembl info
vep <- rbind( # keep only alleles we're examining
  vep[peaks[,.(ID, ALT)], on = .(ID, variant_allele=ALT), nomatch=0],
  vep[peaks[,.(ID, REF)], on = .(ID, variant_allele=REF), nomatch=0]
)
vep <- vep[biotype == "protein_coding"] 
vep <- vep[!(most_severe %in% c("upstream_gene_variant", "downstream_gene_variant"))]
vep <- unique(vep)
vep <- vep[, .(most_severe, consequence_terms, variant_allele, impact, polyphen_prediction, 
               sift_prediction=paste(unique(sift_prediction), collapse=";")), by=.(ID, Gene)]
vep <- unique(vep)

# add to peak information
vep <- vep[, .(ID, Gene, VEP_consequence = most_severe, VEP_impact = impact, polyphen_prediction, sift_prediction)]
vep[sift_prediction == "NA", sift_prediction := NA]
peaks <- merge(peaks, vep, by = "ID", all.x=TRUE)

# For peaks not in genes, find closest coding gene
library(annotables)
genes <- as.data.table(grch37)
genes <- genes[chr %in% 1:22 & biotype == "protein_coding"]
genes[, chr := as.integer(chr)]

to_anno <- peaks[is.na(Gene), .(ID, CHR=`#CHROM`, POS)]
to_anno <- to_anno[order(POS)][order(CHR)]
to_anno <- unique(to_anno)

closest <- foreach(ii = to_anno[,.I], .combine=rbind) %do% {
  var_chr <- to_anno[ii, CHR]
  var_pos <- to_anno[ii, POS]
  closest <- genes[chr == var_chr, .(Gene=symbol, CHR=chr, dist=abs(start - var_pos))][dist == min(dist)]
  closest <- unique(closest)
  cbind(to_anno[ii], Gene = closest$Gene, distance = var_pos - genes[symbol == closest$Gene, start])
}
closest <- closest[abs(distance) < 1e6] # all genes within 1 Mb
peaks[closest, on = .(`#CHROM`=CHR, POS, ID), Gene := i.Gene]

fwrite(peaks, sep="\t", quote=FALSE, file="paper_output/gwas_peak_annotations.txt")

