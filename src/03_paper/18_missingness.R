library(data.table)
library(ukbnmr)
library(foreach)
library(doMC)
library(ggplot2)
library(ggrastr)

# Create output directory
if (!dir.exists("paper_output")) dir.create("paper_output")

# Load data with adjusted values at each step
dat <- fread("data/tech_qc/multistep_adjusted_values.txt")

# Load technical information and filter to samples in UKB raw data
sinfo <- fread("data/tech_qc/sample_information.txt")
sinfo <- sinfo[!(sample_removed) & (in_ukb_raw)]
dat <- dat[sinfo[,.(sample_id, visit)], on = .(sample_id, visit)]

# Compute and plot missigness rates
sample_miss <- dat[,.(N=sum(is.na(adj5)), pct=sum(is.na(adj5))/.N*100),by=.(eid_7439, visit)]
bio_miss <- dat[,.(N=sum(is.na(adj5)), pct=sum(is.na(adj5))/.N*100), by=variable]

# Create density plots
g <- ggplot(sample_miss) + 
  aes(x=pct) +
  geom_density(color="#0000d4") + 
  rasterize(geom_rug(size=0.25, alpha=0.5, color="#d40000"), dpi=1200) +
  scale_x_continuous(name="% biomarkers missing per sample", expand=c(0,0.1)) +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7))
ggsave(width=3.6, height=1.3, file="paper_output/postqc_sample_missingness.pdf")

g <- ggplot(bio_miss) + 
  aes(x=pct) +
  geom_density(color="#0000d4") + 
  rasterize(geom_rug(size=0.25, alpha=0.5, color="#d40000"), dpi=1200) +
  scale_x_continuous(name="% missing per biomarker", expand=c(0,0.1)) +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7))
ggsave(width=3.6, height=1.3, file="paper_output/postqc_biomarker_missingness.pdf")

# Now in raw data
sample_miss <- dat[variable %in% nmr_info[(Nightingale), Biomarker],.(N=sum(is.na(raw)), pct=sum(is.na(raw))/.N*100),by=.(eid_7439, visit)]
bio_miss <- dat[variable %in% nmr_info[(Nightingale), Biomarker],.(N=sum(is.na(raw)), pct=sum(is.na(raw))/.N*100), by=variable]

# Create density plots
g <- ggplot(sample_miss) + 
  aes(x=pct) +
  geom_density(color="#0000d4") + 
  rasterize(geom_rug(size=0.25, alpha=0.5, color="#d40000"), dpi=1200) +
  scale_x_continuous(name="% biomarkers missing per sample", expand=c(0,0.1)) +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7))
ggsave(width=3.6, height=1.3, file="paper_output/raw_sample_missingness.pdf")

g <- ggplot(bio_miss) + 
  aes(x=pct) +
  geom_density(color="#0000d4") + 
  rasterize(geom_rug(size=0.25, alpha=0.5, color="#d40000"), dpi=1200) +
  scale_x_continuous(name="% missing per biomarker", expand=c(0,0.1)) +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=7))
ggsave(width=3.6, height=1.3, file="paper_output/raw_biomarker_missingness.pdf")

