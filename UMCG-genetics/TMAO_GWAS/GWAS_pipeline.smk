##################################################################################################
#####Snakemake pipeline for GWAS analysis of TMAO concentrations in two cohorts LLD and 300OB.####
##################################################################################################



#Steps:
#1. Check/remove family relations
#2. LD pruning, MAF filtering and genotype calling filtering for PCAv
#3. Create PCA plot, VISUAL INSPECTION REQUIRED either include as covariate or don't
#4. Fit lm on Phenotype ~ Age + Sex + (PC) + (kinship) + (Further covariates) + (Fasting) ; Residual transformation
#5. Fit lm on transf(Residuals_4) ~ Genotype ; as an additive model
#6. Outcome formatting



# Important modules
#  module load plink/1.9-foss-2015b
#  Anconda with local isntallation of rvtests 


if config["Filter"] == "yes":
	rule Remove_samples:
		input: 
			remove = config["Remove_list"],
			Genotype_file = 
		output: 
		shell:
			"""
			module load plink/1.9-foss-2015b
			count=$(wc -l {input})
			if [$count > 0]; then
				plink --remove {input.remove}  

			"""
else:
	rule symlink:
		input: Genotype_file = ""
		output:
		shell:
			"""
			ln -s {input.Genotype_file} {output}
			"""
rule filter_PCA:
	input:
		Genotype_file =
	output:
	shell:
		"""
		module load plink/1.9-foss-2015b
		plink --file {input.Genotype_file} --indep-pairwise 50 5 0.3 --maf 0.01 --geno 0.1 --hwe 0.000001 --out {params.prefix}
		"""
rule do_PCA:
	input: rules.oufilter_PCA = 
	output:
	shell:
		"""
		module load plink/1.9-foss-2015b
		plink --file --pca --tabs
		"""
		#Make plot with R?
rule Covariate_model:
	input:
	output:
	shell:
		#R script for lm and transformation
rule SNP_association:
	input:
		Genotypes="",
		Dependent=rules.Covariate_model.output
	output:
	shell:
		"""
		rvtests --inVcf input.vcf --pheno metabo.pheno pheno-name residuals --dosage DS --out output --meta score
		"""
		
