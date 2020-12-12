library(TwoSampleMR)
library(readxl)

##### Exposure #####
## read in the exposure data from a file ##
## Exposure for all of the data in the reference spread sheet ##
Ins_all = read_excel("~/Desktop/STA5934 Statistical Genetics/Final Project/ref_exp_all.xlsx",1)
Ins_all = Ins_all[-which(Ins_all$Chr==6),] #remove the SNPs and proteins encoded by genes with the MHC region (chr6)
Ins_all = format_data(Ins_all, type = "exposure", header = TRUE, phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se", 
                      eaf_col = "MAF", effect_allele_col = "Major_allele", other_allele_col = "Minor_allele", pval_col = "pval")
#rs35233997, rs11292716, rs5784945, rs35400166, rs35365539 are removed due to the lack of required info
Ins_all = clump_data(Ins_all) #clump the data to make sure that the SNPs are conditionally independent using EUR population reference with r2<0.001

##### specificity test #####
specif = matrix(NA, numIV_all, 2)
rownames(specif) = unique(Ins_all$SNP)
colnames(specif) = c("SNP", "NumOfAssProt")
for (i in 1:numIV_all){
  specif[i,1] = UniqueSNP_all[i]
  specif[i,2] = length(which(Ins_all$SNP==UniqueSNP_all[i]))
}
specif = data.frame(specif)
rm = specif$SNP[which(specif$NumOfAssProt>5)] # 3 SNPs associated with more than 5 proteins will be removed
for(i in 1:length(rm)){
  Ins_all = subset(Ins_all, SNP != rm[i])
}
UniqueSNP_all_filter = unique(Ins_all$SNP)
numIV_all_filter = length(unique(Ins_all$SNP)) #226 conditionally independent SNPs remained
numProt_all_filter = length(unique(Ins_all$exposure)) #247 proteins remained 

##### Outcome #####
outcome_final = c("finn-a-F5_SCHZPHR", #Schizophrenia
                  "finn-a-KRA_PSY_SCHIZODEL", #Schizophrenia or delusion
                  "finn-a-F5_SCHIZO", #Schizophrenia, schizotypal and delusional disorders
                  "finn-a-F5_SCHIZOTYP" #Schizotypal disorder
)
out_finn = extract_outcome_data(snps = Ins_all$SNP, outcomes = outcome_final) #proxies with LD(r2>0.8)
#184 snps for each outcome id are extracted when extract all together

##### Harmonise the data #####
## for wald ratio and ivw ##
scz_mrdat_all = harmonise_data(exposure_dat = Ins_all, outcome_dat = out_finn, action = 2) #harmonise data

## for mr-egger ##
flip = which(scz_mrdat_all$beta.exposure < 0)
scz_mrflip_all = scz_mrdat_all
scz_mrflip_all$beta.exposure[flip] = -1 * scz_mrflip_all$beta.exposure[flip]
scz_mrflip_all$beta.outcome[flip] = -1 * scz_mrflip_all$beta.outcome[flip]

##### MR analyses #####
## wald ratio test for exposures associated with single snps sequentially ##
mr_seq_wald = mr(scz_mrdat_all, method_list=c("mr_wald_ratio"))
num_single = length(unique(mr_seq_wald$exposure)) #191 proteins associated with single SNP
pthr_wald = 0.05/num_single
mr_seq_wald_causal = mr_seq_wald[which(mr_seq_wald$pval < pthr_wald),] # no causal phynotypes identified 

## IVW (fixed effects) test for exposures associated with multiple snps sequentially ##
mr_seq_ivw = mr(scz_mrdat_all, method_list=c("mr_ivw_fe"))
num_mul = length(unique(mr_seq_ivw_mre$exposure)) #13 proteins associated with 2 SNPs
pthr_ivw = 0.05/(num_mul*2)
length(unique(mr_seq_ivw$exposure))
mr_seq_ivw_causal = mr_seq_ivw[which(mr_seq_ivw$pval < pthr_ivw),] # no causal phynotypes identified 

## IVW (multiplicative random effects) test for exposures associated with multiple snps sequentially
mr_seq_ivw_mre = mr(scz_mrdat_all, method_list=c("mr_ivw_mre"))
mr_seq_ivw_mre_causal = mr_seq_ivw_mre[which(mr_seq_ivw_mre$pval < pthr_ivw),] 

## MR-Egger ##
# can not run mr-egger and the pleiotropy test either, the reason might be there's not enough SNPs for this mehtod
# according to my test, mr-egger works for at lease three SNPs, our phynotypes and traits are associated with at most 2 SNPs
# mr_seq_egger = mr(scz_mrflip_all, method_list=c("mr_egger_regression")) # no result

##### Sensitivity analyses #####
## heterogeneity test wtih multiplicative random effects IVW ##
mr_hetero <- mr_heterogeneity(scz_mrdat_all, method_list=c("mr_ivw") ) 
subset(mr_hetero, Q_pval < 0.1)

## Bidirectional MR ##
#harmonise again
bi_exp_scz = out_finn
bi_exp_scz = bi_exp_scz[,-c(13,14,16,17,18,19,20,21,22,23)]
colnames(bi_exp_scz) = c("SNP", "chr", "pos", "beta.exposure", "se.exposure", "samplesize.exposure", "pval.exposure", "eaf.exposure", 
                         "effect_allele.exposure", "other_allele.exposure", "exposure", "id.exposure", "mr_keep.exposure")
bi_out_prot = Ins_all
bi_out_prot = bi_out_prot[,-10]
colnames(bi_out_prot) = c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", 
                          "outcome", "mr_keep.outcome", "id.outcome")
bi_dir_dat = harmonise_data(exposure_dat = bi_exp_scz, outcome_dat = bi_out_prot, action = 2)

bimr_seq_wald = mr(bi_dir_dat, method_list=c("mr_wald_ratio"))
bimr_seq_wald_causal = bimr_seq_wald[which(bimr_seq_wald$pval < pthr_wald),]
#dim(bimr_seq_wald)[1] == dim(bimr_seq_wald_causal)[1] is true, meaning all proteins showed reverse causalitiy 

bimr_seq_ivw_mre = mr(bi_dir_dat, method_list=c("mr_ivw_mre"))
bimr_seq_ivw_mre_causal = bimr_seq_ivw_mre[which(bimr_seq_ivw_mre$pval < pthr_ivw_mre),]