library(data.table)
library(ukbnmr)
library(MASS)
library(readstata13)
library(ggplot2)
source("src/utilities/logit.R")

# Create output directory
if (!dir.exists("diagnostic_plots")) dir.create("diagnostic_plots")

# Adjust post-qc data for age, sex, and bmi, applied to all biomarkers
# and also with biomarker rederivation so we can compare rederivation
# to no rederivation in a realistic scenario where we expect much larger
# changes after adjustment than those from removing technical variation.
dat <- fread("data/tech_qc/postqc_excluding_outlier_plates.txt")

# For blind duplicates, filter to samples kept by UK Biobank
sinfo <- fread("data/tech_qc/sample_information.txt")
dat <- dat[sinfo[(in_ukb_raw) & !(sample_removed), .(sample_id, visit)], on = .(sample_id, visit)]

# Load in information about participant age and BMI from departmental
# curated dataset.
ext_rep <- read.dta13("data/ceu_curated_phenotypes/20210901/STATA/repeats.dta")
setDT(ext_rep)
pheno <- ext_rep[,.(eid_7439=idno, visit=repno, age=ages_rep, bmi=bmi_rep)]
pheno <- pheno[visit %in% c(0,1)]
pheno[, visit := ifelse(visit == 0, "Main Phase", "Repeat Assessment")]

# Load in participant sex
sex <- read.dta13("data/ceu_curated_phenotypes/20210901/STATA/analysis.dta")
setDT(sex)
sex <- sex[,.(eid_7439=idno, sex)]

# Add to dat
dat[, eid_7439 := as.character(eid_7439)]
dat <- dat[pheno, on = .(eid_7439, visit), nomatch = 0]
dat <- dat[sex, on = .(eid_7439), nomatch = 0]

# Melt to long format
dat <- melt(dat, id.vars=c("eid_7439", "sample_id", "visit", "age", "sex", "bmi"), measure.vars=ukbnmr::nmr_info[(Nightingale), Biomarker], value.name="postqc")

# Percentages will be logit transformed instead of log transfomred
dat <- dat[!is.na(postqc)]
dat[variable %in% ukbnmr::nmr_info[Units == "%", Biomarker], postqc := postqc/100]

# Log/logit offset for 0 values
lower_lim <- dat[!is.na(postqc),.(min=min(postqc), non_zero=min(postqc[postqc != 0])),by=variable]
lower_lim[, log_offset := ifelse(min == 0, non_zero / 2, 0)]
dat[lower_lim, on = .(variable), postqc := postqc + log_offset]

# Log / logit transformation
dat[variable %in% ukbnmr::nmr_info[Units == "%", Biomarker], log_postqc := logit(postqc)]
dat[variable %in% ukbnmr::nmr_info[Units != "%", Biomarker], log_postqc := log(postqc)]

# Function for converting a set of numbers to a factor, ordering levels by group size
factor_by_size <- function(x) {
  group_sizes <- as.data.table(table(x))
  group_sizes <- group_sizes[order(-N)]
  factor(x, levels=group_sizes$x)
}

# Adjust for age, sex, and bmi:
dat <- dat[!is.na(sex) & !is.na(age) & !is.na(bmi)]
dat[, adj := rlm(postqc ~ age + factor_by_size(sex) + log(bmi))$residuals, by=variable]

# Rescale the adjusted residuals back to absolute units
rescale <- function(log_adj, log_raw) {
  # The residuals from any regression are defined as the difference between the
  # observed independent variable (e.g. biomarker concentrations) and the parameter
  # estimated by the regression. A resulting key property is that the distribution
  # of the residuals is centred on 0 for this estimated parameter. In the case of
  # robust linear regression, this parameter is an estimate of the mean that is
  # robust to outliers. As a consequence of the way residuals are defined, their
  # distribution can be scaled to match the distribution of the independent variable
  # by giving it the same estimated mean. For robust linear regression, residuals
  # can be returned to the same scale as the observed independent variable by
  # estimating the mean of the observed independent variable using robust linear
  # regression and adding it to the residuals.
  adj <- log_adj + as.vector(coef(rlm(log_raw ~ 1)))[1]
  return(adj)
}
dat[, adj := rescale(adj, log_postqc), by=.(variable)]
dat[variable %in% ukbnmr::nmr_info[Units == "%", Biomarker], adj := invlogit(adj)]
dat[variable %in% ukbnmr::nmr_info[Units != "%", Biomarker], adj := exp(adj)]
dat[lower_lim, on = .(variable), adj := adj - log_offset]

# Add offset so that there are no negative values
shift <- dat[,.(right_shift=-pmin(0, min(adj))), by=variable]
lower_lim <- lower_lim[shift, on = .(variable)]
dat[lower_lim, on = .(variable),  adj := adj + right_shift]

# Percentages back to 0-100 instead of 0 to 1
dat[variable %in% ukbnmr::nmr_info[Units == "%", Biomarker], adj := adj * 100]

# We want to create two sets for comparison - one with composite and derived
# biomarkers computed from their adjusted parts, and one where the biomarkers 
# have been adjusted through regression.
non_deriv <- dat[variable %in% ukbnmr::nmr_info[Type == "Non-derived", Biomarker]]
non_deriv <- dcast(non_deriv, eid_7439 + sample_id + visit ~ variable, value.var="adj")
non_deriv <- recompute_derived_biomarkers(non_deriv)
non_deriv <- melt(non_deriv, id.vars=c("eid_7439", "sample_id", "visit"), value.name="adj")

dat <- dat[, .(eid_7439, sample_id, visit, age, sex, bmi, variable, adj_no_rederive=adj)]
dat[non_deriv, on = .(eid_7439, sample_id, visit, variable), adj := i.adj]

# Write out
if (!dir.exists("data/age_sex_bmi_adj")) dir.create("data/age_sex_bmi_adj")
fwrite(dat, sep="\t", quote=FALSE, file="data/age_sex_bmi_adj/adjusted_vs_recomputed.txt")

# Show difference between methods:
if (!dir.exists("diagnostic_plots")) dir.create("diagnostic_plots")

gg_dt <- dat[variable %in% ukbnmr::nmr_info[Type != "Non-derived" & (Nightingale), Biomarker]]
gg_dt <- gg_dt[!is.na(adj) & !is.na(adj_no_rederive)]
gg_dt <- gg_dt[order(tolower(variable))]
gg_dt[, variable := factor(variable, levels=unique(variable))]

g <- ggplot(gg_dt) +
  aes(x=adj, y=adj_no_rederive) +
  geom_point(shape=21, color="black", fill="white", size=0.45, stroke=0.2) +
  geom_abline(intercept=0, slope=1, linetype=2, color="red", size=0.4) +
  scale_x_continuous(name="Derived and composite biomarkers computed from adjusted parts") +
  scale_y_continuous(name="Biomarker adjusted for technical covariates, age, sex, and BMI") +
  facet_wrap(~ variable, ncol=12, scales="free") +
  theme_bw() +
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=8)
  )
ggsave(g, width=16.8, height=16.8, units="in", file="diagnostic_plots/adjusted_vs_from_parts_age_sex_bmi.png")

