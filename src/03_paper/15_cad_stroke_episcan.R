library(data.table)
library(ukbnmr)
library(survival)
library(readstata13)
library(foreach)
library(ggplot2)
library(ggrastr)
library(cowplot)

options(ggrastr.default.dpi=1200)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Load technical information and filter to samples in UKB raw data
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(sample_id, visit)], on = .(sample_id, visit)]

# Filter to participants at baseline assessment
dat <- dat[visit == "Main Phase"]

# Load in curated CAD and Stroke endpoints
cad <- fread("data/my_curated_phenotypes/endpoints/output/Hard_CAD/Hard_CAD_events_and_followup.txt")
cad <- cad[visit == 0]
dat[cad, on = .(eid_7439=eid), c("CAD_event", "CAD_follow", "CAD_prevalent") := .(inci_event, inci_followup, prev_event)]

stroke <- fread("data/my_curated_phenotypes/endpoints/output/CEU_Stroke/CEU_Stroke_events_and_followup.txt")
stroke <- stroke[visit == 0]
dat[stroke, on = .(eid_7439=eid), c("Stroke_event", "Stroke_follow", "Stroke_prevalent") := .(inci_event, inci_followup, prev_event)]

# Log or logit transform biomarkers
source("src/utilities/logit.R")

dat[!is.na(raw) & !is.finite(raw), raw := NA]
dat[variable %in% nmr_info[Units == "%", Biomarker], raw := raw/100]
lim_raw <- dat[!is.na(raw), .(min=min(raw), non_zero=min(raw[raw != 0]), non_one=max(raw[raw != 1]), max=max(raw)), by=variable]
lim_raw[variable %in% nmr_info[Units != "%", Biomarker], offset := ifelse(min == 0, non_zero / 2, 0)]
lim_raw[variable %in% nmr_info[Units == "%", Biomarker], offset := fcase(
  min == 0 & max == 1, 0, # cannot offset
  min == 0, pmin(non_zero / 2, (1 - max)/2), # make sure offset doesnt make % > 100
  max == 1, -pmin(min / 2, (1 - non_one)/2), # make sure offset doesnt make % < 0
  default = 0)]
dat[lim_raw, on = .(variable), raw := raw + offset]
dat[variable %in% nmr_info[Units == "%", Biomarker], raw := logit(raw)]
dat[variable %in% nmr_info[Units != "%", Biomarker], raw := log(raw)]

dat[!is.na(adj5) & !is.finite(adj5), adj5 := NA]
dat[variable %in% nmr_info[Units == "%", Biomarker], adj5 := adj5/100]
lim_adj5 <- dat[!is.na(adj5), .(min=min(adj5), non_zero=min(adj5[adj5 != 0]), non_one=max(adj5[adj5 != 1]), max=max(adj5)), by=variable]
lim_adj5[variable %in% nmr_info[Units != "%", Biomarker], offset := ifelse(min == 0, non_zero / 2, 0)]
lim_adj5[variable %in% nmr_info[Units == "%", Biomarker], offset := fcase(
  min == 0 & max == 1, 0, # cannot offset
  min == 0, pmin(non_zero / 2, (1 - max)/2), # make sure offset doesnt make % > 100
  max == 1, -pmin(min / 2, (1 - non_one)/2), # make sure offset doesnt make % < 0
  default = 0)]
dat[lim_adj5, on = .(variable), adj5 := adj5 + offset]
dat[variable %in% nmr_info[Units == "%", Biomarker], adj5 := logit(adj5)]
dat[variable %in% nmr_info[Units != "%", Biomarker], adj5 := log(adj5)]

# Load in information about additional covariates from departmental curated dataset
ext <- read.dta13("data/ceu_curated_phenotypes/20210302/STATA/analysis.dta")
setDT(ext)
ext <- ext[,.(eid_7439=idno, age=ages, sex, lipid_med = lipdbin)]
dat[, eid_7439 := as.character(eid_7439)]
dat <- dat[ext, on = .(eid_7439), nomatch=0]

# Remove participants taking lipid lowering medication
dat <- dat[lipid_med != "Current"]

# Get hazard ratios
hrs <- rbind(idcol="endpoint",
  CAD = foreach(type = c("raw", "adj5"), .combine=rbind) %:% foreach(var = unique(dat$variable), .combine=rbind) %dopar% {
    this_dat <- dat[variable == var & CAD_prevalent == 0]
    setnames(this_dat, type, "value")
    if (this_dat[!is.na(value), .N] == 0) return(NULL)

    cx <- coxph(Surv(CAD_follow, CAD_event) ~ scale(value) + age + factor(sex), data=this_dat)
    ci <- confint(cx)
    cf <- coef(summary(cx))
  
    data.table(samples=this_dat[,.N], cases=this_dat[,sum(CAD_event)], variable = var, dataset = type, 
               logHR = cf[1,1], SE = cf[1,3], HR=cf[1,2], L95=exp(ci[1,1]), U95=exp(ci[1,2]), P=cf[1,5])
  },
  Stroke = foreach(type = c("raw", "adj5"), .combine=rbind) %:% foreach(var = unique(dat$variable), .combine=rbind) %dopar% {
    this_dat <- dat[variable == var & Stroke_prevalent == 0]
    setnames(this_dat, type, "value")
    if (this_dat[!is.na(value), .N] == 0) return(NULL)

    cx <- coxph(Surv(Stroke_follow, Stroke_event) ~ scale(value) + age + factor(sex), data=this_dat)
    ci <- confint(cx)
    cf <- coef(summary(cx))
  
    data.table(samples=this_dat[,.N], cases=this_dat[,sum(Stroke_event)], variable = var, dataset = type,
               logHR = cf[1,1], SE = cf[1,3], HR=cf[1,2], L95=exp(ci[1,1]), U95=exp(ci[1,2]), P=cf[1,5])
  }
)

# Compare P-value distributions
gg_dt <- dcast(hrs, endpoint + variable ~ dataset, value.var="P")
gg_dt <- gg_dt[!is.na(raw)] # drops extended biomarkers not in raw data
gg_dt[, raw := -log10(raw)]
gg_dt[, adj5 := -log10(adj5)]

g <- ggplot(gg_dt) + 
  aes(x=raw, y=adj5) +
  rasterize(geom_point(shape=21, color="white", fill="#DE6078", stroke=0.1, size=1)) +
  geom_abline(intercept=0, slope=1, linetype=2, color="black", size=0.3) +
  xlab("Raw -log10 P-values") + 
  ylab("Post-QC -log10 P-values") +
  facet_wrap(~ endpoint, scales="free") +
  theme_bw() + 
  theme(axis.text.x=element_text(size=6), axis.title.x=element_text(size=7),
        axis.text.y=element_text(size=6), axis.title.y=element_text(size=7),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(), strip.text=element_text(size=7)
  )

ggsave(g, width=3.6, height=2, file="paper_output/epi_compare.pdf")

# Output hazard ratios
hrs <- dcast(hrs, variable ~ dataset + endpoint, value.var=c("HR", "L95", "U95", "P"))
hrs <- hrs[order(P_adj5_CAD), .(variable,
  HR_adj5_CAD, L95_adj5_CAD, U95_adj5_CAD, P_adj5_CAD,
  HR_raw_CAD, L95_raw_CAD, U95_raw_CAD, P_raw_CAD,
  HR_adj5_Stroke, L95_adj5_Stroke, U95_adj5_Stroke, P_adj5_Stroke,
  HR_raw_Stroke, L95_raw_Stroke, U95_raw_Stroke, P_raw_Stroke
)]
fwrite(hrs, sep="\t", quote=FALSE, file="paper_output/epi_scan.txt")
