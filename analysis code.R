# MR analysis
library(ieugwasr)
library(TwoSampleMR)
library(data.table)
library(plyr)
library(dplyr)
library(gwasglue)
library(TwoSampleMR)
library(vroom)
data<-vroom::vroom("\\~")
exposure<-subset(data, pval.exposure<"Threshold for P-value")

data<-read_exposure_data(filename="\\~",
                         sep = "",
                         snp_col = "",
                         beta_col = "",
                         se_col = "",
                         pval_col = "",
                         effect_allele_col="",
                         other_allele_col = "",
                         eaf_col = "",
                         phenotype_col = "",
                         id_col = "",
                         samplesize_col = "",
                         chr_col="", 
                         pos_col = "",
                         clump=F)
# exposure_dat<- clump_data(exposure_dat, clump_kb=10000, clump_r2=0.001,clump_p1 = 1,clump_p2 = 1,pop = "EUR")
out<- ld_clump(dplyr::tibble(rsid=exposure_dat$SNP,pval=exposure_dat$pval.exposure,id=exposure_dat$id.exposure),
               clump_kb = 10000,
               clump_p = 1,
               clump_r2 = 0.001,
               plink_bin = plinkbinr::get_plink_exe(),
               bfile = "\\~",
               pop = "EUR")
keep<-paste(exposure_dat$SNP,exposure_dat$id.exposure)%in%paste(out$rsid,out$id)
exposure_dat<-exposure_dat[keep,]

exposure_dat$R2<-(2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)/(2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)+2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)))
exposure_dat$F<-exposure_dat$R2*(exposure_dat$samplesize-2)/(1-exposure_dat$R2)
outTab<-exposure_dat[as.numeric(exposure_dat$F)>10,]

exposure_dat<-read.csv("\\~.csv")
outcomedata=read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "\\~",
  sep = "",
  snp_col = "",
  beta_col = "",
  se_col = "",
  effect_allele_col = "",
  other_allele_col = "",
  eaf_col = "",
  pval_col = "",
  ncase_col = "",
  ncontrol_col = "",
  samplesize_col = "",
  id_col = "",
  chr_col = "",
  pos_col = "")
dat=harmonise_data(exposure_dat, outcomedata)
mrResult=mr(dat,
            method_list = c("mr_ivw",
                            # "mr_ivw_fe",
                            "mr_egger_regression",
                            "mr_weighted_median",
                            "mr_simple_mode",
                            "mr_weighted_mode"))
mrTab=generate_odds_ratios(mrResult)
heterTab=mr_heterogeneity(dat)
pleioTab=mr_pleiotropy_test(dat)
presso=run_mr_presso(dat)


# Meta analysis
library(meta)
dat<-read.csv("\\~.csv")
# mata_Result<-meta::metagen(log(dat$or),(log(dat$or_uci95)-log(dat$or_lci95))/(2*1.96),sm="OR",data=dat,studlab=dat$id.outcome)
mata_Result<-meta::metagen(log(dat$or),(log(dat$or_uci95)-log(dat$or_lci95))/(2*1.96),sm="OR",data=dat,studlab=dat$id.outcome, control = list(maxiter=10000,stepadj=0.5))
meta_R<-summary(mata_Result)
Egger<-metabias(mata_Result, k.min = 5, method.bias = "Egger")
Begg<-metabias(mata_Result, k.min = 5, method.bias = "Begg", correct =T, plotit = F)


# LDSC analysis
library(ldscr)
rg_res <- ldsc_rg(
  munged_sumstats = list("A" = "\\~",
                         "B" = "\\~"),
  ancestry = "EUR"
)