# S. Andreu-Sanchez, 7-May-2024

#Prior to running script: ml RPlus/4.2.1-foss-2022a-v22.12


.libPaths(c("/groups/umcg-llnext/tmp01/umcg-sandreusanchez/HLA_imputation/scripts/Package", .libPaths()))
library(tidyverse)
library(HIBAG)

#Sergio's comments:
#HIBAG is a software that allows imputation of HLA class I (A,B,C) and class II (HLA-DRB1, HLA-DQB1, HLA-DPB1, HLA-DQA1) from genotype data.
#Genotype data should be stored in Plink 1 format (Bed, bim, fam) files to use the script. Ideally these files would QC-ed, not imputed, genotypes. Since the algorithm takes observed genotypes as truth, and does not use genotypes likelihoods from input genotypes (allowing to take into account uncertainty introduced by genome-wide imptuation). In practice, for the chip that I tested (the LL-next data), there were no genotype SNPs that matched the SNPs that HIBAG uses, so I just used the QC-ed genome-wide imputed genotypes. 
#HIBAG is already pre-trained in European data (to train again, MHC needs to be known, and steps can be followed from the HIBAG documentation). The algorithm is an ensemble method, so there are a bunch of different models (100) fitted in different subsets of training data with slighly different SNPs used for prediction. Each of the predictors is fitted with an EM algortihm, outputting posterior probabilites for each possible haplotype. The results are the average of the posteriors for each of the data subsets.  




##Pre-trained model, downloaded from https://hibag.s3.amazonaws.com/download/hlares_param/European-HLA4.html on 7/May/2024. "wget https://hibag.s3.amazonaws.com/download/hlares_param/European-HLA4.html"
model.list <- get(load("/groups/umcg-llnext/tmp01/umcg-sandreusanchez/HLA_imputation/European-HLA4-hg19.RData"))

##Get your files. This might take several minutes 

#Old imputed data
BED = "/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.info_filtered.bed"
BIM = "/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.info_filtered.bim"
FAM = "/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/all_chr.info_filtered.fam"
#New unimputed data (no SNPs in the region)
#BED = "/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotype_QC/output/8_second_iteration/final_filtered_genotypes/llnext_batch12.bed"
#BIM = "/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotype_QC/output/8_second_iteration/final_filtered_genotypes/llnext_batch12.bim"
#FAM = "/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotype_QC/output/8_second_iteration/final_filtered_genotypes/llnext_batch12.fam"


yourgeno <- hlaBED2Geno(bed.fn=BED, fam.fn=FAM, bim.fn=BIM, import.chr="6", assembly="hg19")


Do_imputation = function(model.list, Allele = "A" ){
	if (!Allele %in% Alleles){ print( paste0("Allele not in allele list: ", pate(Alleles, collapse=",")    ) ) }
	model <- hlaModelFromObj(model.list[[Allele]])

	#Check number of matching positions
	N_SNP_model = model$snp.position
	SNP_match = model$snp.position[model$snp.position  %in%  yourgeno$snp.position]

	#Do prediction
	pred.guess <- hlaPredict(model, yourgeno, type="response+prob")


	#Each sample is in a row. Only one entry per sample. The sample with the highest posterior is chosen. The probability given is the posterior. 
	#Matching proportion is a measure or proportion describing how the SNP profile matches the SNP haplotypes observed in the training set, i.e., the likelihood of SNP profile in a random-mating population consisting of training haplotypes. It is not directly related to confidence score, but a very low value of matching indicates that it is underrepresented in the training set.
	pred.guess$value %>% as_tibble() %>% mutate(Allele = Allele)  -> ImputedGenotypes
	#Each sample is a column. Each allele is a row. The sum of posteriors per sample sum up to 1.
	pred.guess$postprob %>% as.data.frame() %>% rownames_to_column("Allele") %>% as_tibble() %>% mutate(Allele = Allele)  -> Posteriors

	return(list(ImputedGenotypes, Posteriors))
}


print("Fitting models for each possible Allele")

Alleles = c("A",    "B",    "C",    "DRB1", "DQA1", "DQB1", "DPB1")



Posteriors = tibble()
Haplotypes = tibble()
for (Alle in Alleles){
	print( paste0("Iterating allele: ", Alle) )
	Result = Do_imputation(model.list, Alle)

	Posteriors %>% rbind(Result[[2]]) -> Posteriors
	Haplotypes %>% rbind(Result[[1]]) ->  Haplotypes

}
Haplotypes %>% unite(Haplotype, allele1, allele2, sep = ";") %>% select(-c(matching, prob)) %>% spread(Allele, Haplotype) -> Haplotypes2


Output_haplo =  "/groups/umcg-llnext/tmp01/umcg-sandreusanchez/HLA_imputation/LLNEXT_imputedv1_HLA_HaplotypeImputationHIBAG.tsv"
Output_haplo_wide =  "/groups/umcg-llnext/tmp01/umcg-sandreusanchez/HLA_imputation/LLNEXT_imputedv1_HLA_HaplotypeImputationHIBAG_final.tsv"
Output_posterior =  "/groups/umcg-llnext/tmp01/umcg-sandreusanchez/HLA_imputation/LLNEXT_imputedv1_HLA_HaplotypeImputationHIBAG_posteriors.tsv"

write_tsv(Haplotypes, Output_haplo)
write_tsv(Posteriors, Output_posterior)
write_tsv(Haplotypes2, Output_haplo_wide)



