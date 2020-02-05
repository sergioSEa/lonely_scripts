from pathlib import Path
configfile: "config.yaml"

Welcome_message = """
snakemake --latency-wait 60 --rerun-incomplete  --jobs 99 --keep-going --cluster  'sbatch -A PROJECT  -t {cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.c} --error={cluster.error} --job-name={cluster.name} --output={cluster.output}'  --cluster-config cluster.json --snakefile Snake_GWAS.smk  --cluster-status 'python slurm-status.py'
"""
print(Welcome_message)

localrules: Choose_phenotypes, Split_Male_Female

Input_path = config["Input"]
#Cohort = config["Cohort"]

#Input_path="/groups/umcg-lld/tmp03/umcg-akurilshchikov/LLD_HRC"

Outputs = []

Outputs.append("Output/Plot_PCA.png")
Outputs.append("Output/plink.genome")
#Final output to be appended

Inputs = []
Intermediates = []
List_names = []
for F in Path(Input_path).glob("*.vcf.gz"):
	Inputs.append(F)
	Inter =  "Filtered/" + F.stem.split(".")[0] + ".vcf"
	List_names.append(F.stem.split(".")[0])
	Intermediates.append(Inter)

	
rule all:
	input: 
		Outputs,
		expand("Output/GWAS_results_{Pheno}{Gender}.MetaScore.assoc.gz", Pheno=["TMAO","Choline","Betaine","L-Carnitine"],Gender=["1","2","3"]),
		expand("{Pheno}.{Study}.{Gender}.done",  Pheno=["TMAO","Choline","Betaine","L-Carnitine"],Gender=["1","2","3"], Study=config["Cohort"])
		
rule Filter_individuals:
	input: 
		vcf = Input_path + "/{Name}.dose.vcf.gz",
		remove = config["Remove_individuals"]
	output: "Dataset/{Name}.vcf"
	threads: 6
	shell:
		"module load plink/1.9-foss-2015b ;\n"
		"module load BCFtools/1.7-foss-2015b ;\n"
		"bcftools convert  --samples-file {input.remove} --output {output} --output-type v --threads {threads} {input.vcf}"

rule Filters_for_PCA:
	input: 	rules.Filter_individuals.output
	output:	"Filtered/{Name}.vcf"
	threads: 6
	shell:
		"module load plink/1.9-foss-2015b ;\n"
		"plink --vcf {input} --maf 0.01 --geno 0.1 --hwe 0.000001 --out Filtered/{wildcards.Name} --recode vcf"

rule merge_files:
	input: Intermediates
	output: 
		M ="Filtered/Merged.vcf",
		F = touch("merged_Completed")
	shell:
		"module load BCFtools/1.7-foss-2015b ; \n"
		"python3 scripts/Merge.py Filtered {output.M}"
rule LD_Pruning:
	input: 
		flag = rules.merge_files.output.F,
		i = rules.merge_files.output.M
	output: "Filtered/Merged_pruned.bed"
	shell:
		"module load plink/1.9-foss-2015b ;\n"
		"plink --vcf {input.i}  --indep-pairwise 1000 50 0.3 --out PCA_pruned ; \n"
		"plink --vcf {input.i} --extract PCA_pruned.prune.in --make-bed --out Filtered/Merged_pruned"

rule PCA:
	input: rules.LD_Pruning.output
	output: "Output/Plot_PCA.png"
	shell:
		"module load plink/1.9-foss-2015b ;\n"
		"module load Python/2.7.11-foss-2015b ;\n"
		"module load matplotlib/1.5.1-foss-2015b-Python-2.7.11 ;\n"
		"plink --bfile Filtered/Merged_pruned --pca ;\n"
		"python scripts/Plo_pca.py {output}"

rule remove_imputed:
	input:  expand("Dataset/{Name}.vcf", zip, Name= List_names)
	output:	l = "Genotyped_list.txt",
 		final = "Genotyped.vcf"
	shell:
		"module load plink/1.9-foss-2015b ;\n"
		"touch  Genotyped_list.txt ; \n"
		"array=( {input} ) ; \n"
		'for i in "${{array[@]}}" ; do \n'
			'cat "${{i}}" | grep "PASS;GENOTYPED" | cut -f3 >> Genotyped_list.txt ;\n'
		"done ;\n"
		"plink --vcf Filtered/Merged.vcf --recode vcf --out Genotyped --extract  Genotyped_list.txt"

rule IBD:
	input: rules.remove_imputed.output.final
	output: "Output/plink.genome"
	shell:
		"module load plink/1.9-foss-2015b ;\n"
		"plink --vcf Genotyped.vcf --genome --out Output/plink ;\n"
		
rule find_related:
	input: rules.IBD.output
	output:
		related = "remove_related.txt",
		List = "Output/related.tsv"

	shell:
		"python scripts/Check_IB.py > {output.related} ;\n"
		
rule Filter2:
	input:
		Remove_IBD = rules.find_related.output.related,
		vcf = "Dataset/{Name}.vcf"
	output:
		Remove = temp("Remove2_{Name}.txt"), 
		F = "Filter_GWAS/{Name}.vcf"
	threads: 6
	shell:
		"module load BCFtools/1.7-foss-2015b ;\n"
		"cat {input.Remove_IBD} remove_PCA_outliers.txt | sort | uniq > {output.Remove} ;\n"
		"bcftools convert  --samples-file ^{output.Remove} --output {output.F}  --output-type v --threads {threads} {input.vcf}"

rule Merge2:
	input: expand("Filter_GWAS/{Name}.vcf", zip, Name= List_names)
	output: "Filter_GWAS/Merged.vcf"
	shell:
		"module load BCFtools/1.7-foss-2015b ; \n"
		"module load Python/3.6.3-foss-2015b ; \n"
		"python scripts/Merge.py Filter_GWAS {output}"
rule Prepare_covariates:
	input: 
		PCA = rules.PCA.output,
		Covariates = config["Covariates"]
	output: "Transformed_covariates.tsv"
	shell:	
		"Rscript scripts/Transformation.R {input.Covariates}" 
rule Choose_phenotype:
	input: rules.Prepare_covariates.output
	output:	
		Phenofile = "Output/Phenotype.phe",
		Covariates = "Output/Covariates.phe",
	shell:
		"python Prepare_phenotype_input.py"
rule Split_Male_Female:
	input: 
		Phenotype = rules.Choose_phenotype.output.Phenofile,
		Covariates = rules.Choose_phenotype.output.Covariates
	output:
		Phenofile = "Output/Phenotype{Gender}.phe",
		Covariates= "Output/Covariates{Gender}.phe"
	run:
		if wildcards.Gender != "3": 
			shell("""
			awk -F" " '$5 == "{wildcards.Gender}" || $5 == "sex" {{ print  }}' {input.Phenotype} > {output.Phenofile}
			awk -F" " '$5 == "{wildcards.Gender}" || $5 == "sex" {{ print  }}' {input.Covariates} > {output.Covariates}
			""")
		else:
			shell("""
			cp {input.Phenotype} {output.Phenofile}
			cp {input.Covariates} {output.Covariates}
			""")

		
rule Test:
	input: 
		Genotype = rules.Merge2.output,
		Phenotype =  "Output/Phenotype{Gender}.phe",
		Covariates =  "Output/Covariates{Gender}.phe"
	output: 
		O = "Output/GWAS_results_{Pheno}{Gender}.MetaScore.assoc.gz",
		t = "Output/GWAS_results_{Pheno}{Gender}.MetaScore.assoc.gz.tbi"
	params: 
		rvtests = "rvtest/executable/rvtest"
	threads: 6
	run:
		if wildcards.Gender == "3":
			shell("{params.rvtests} --inVcf {input.Genotype} --pheno {input.Phenotype} --pheno-name {wildcards.Pheno} --dosage DS --out Output/GWAS_results_{wildcards.Pheno}{wildcards.Gender} --covar {input.Covariates} --covar-name bmi,sex,age --meta score --numThread 6")
		else:
			shell("{params.rvtests} --inVcf {input.Genotype} --pheno {input.Phenotype} --pheno-name {wildcards.Pheno} --dosage DS --out Output/GWAS_results_{wildcards.Pheno}{wildcards.Gender} --covar {input.Covariates} --covar-name bmi,age --meta score --numThread 6")

rule Get_r2:
	input: rules.Merge2.output
	output: 
		R2 = "Filter_GWAS/R2_values.txt",
		T = temp("Filter_GWAS/temp")
	shell:
		"""
		grep -v "#" {input}  | cut -f 8 > {output.T}
		sed 's/.*R2=//' {output.T} > {output.R2}
		"""
rule Final_format:
	input: 
		Results = rules.Test.output.O,
		TBI = rules.Test.output.t,
		R2 = rules.Get_r2.output.R2
	output:
		End_flag = touch("{Pheno}.{Study}.{Gender}.done")
	run:
		import datetime
		Metabo_dic = {"Betaine":"BETAINE", "Choline":"CHOLINE", "L-Carnitine":"CARNITINE", "TMAO":"TMAO"}
		Gender_dic = {"1":"MALES","2":"FEMALES","3":"ALL"}

		dt = datetime.datetime.today()
		DT = str(dt.year) + str(dt.month) + str(dt.day)

		Output_name = "Output/{Metabo}.HRC.{Study}.EA.{Gender}.{DATE}.txt".format(Metabo=Metabo_dic[wildcards.Pheno],Study=wildcards.Study,Gender=Gender_dic[wildcards.Gender],DATE=DT)

		shell("bash scripts/Process_output.sh {input.Results} {input.R2} {Output_name} {wildcards.Pheno}_{wildcards.Gender}")
		

#arule Final_format:
#	input:
#		rules.Test.output
#	output: touch("{Pheno}.{Study}.{Gender}.done")
#	shell:
#		"module load Python/3.6.3-foss-2015b ;\n"
#		"python scripts/Format_output2.py {input} {wildcards.Pheno} {wildcards.Gender} {wildcards.Study}"
