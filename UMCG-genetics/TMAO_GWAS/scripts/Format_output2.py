# -*- coding: utf-8 -*-
import sys
import gzip
import datetime
from subprocess import Popen, PIPE

Input_file = sys.argv[1]
Metabo = sys.argv[2]
Gender = sys.argv[3]
Study = sys.argv[4]

#["TMAO","Choline","Betaine","y-butyrobetaine","L-Carnitine"
Metabo_dic = {"Betaine":"BETAINE", "Choline":"CHOLINE", "L-Carnitine":"CARNITINE", "TMAO":"TMAO"}
Gender_dic = {"1":"MALES","2":"FEMALES","3":"ALL"}

dt = datetime.datetime.today()
DT = str(dt.year) + str(dt.month) + str(dt.day)


def Check_vcf(Chr,SNP):
	File = "Filter_GWAS/chr{X}.vcf".format(X=Chr)
	command = """ awk '$2 == "{P}" {{print $8 ; exit}}' {F}""".format(P=SNP, F=File)
	O= Popen(command, shell=True,stdout=PIPE, stderr=PIPE)
	Out = O.communicate()[0].decode('UTF-8')
	R = Out.split("R2=")[1].rstrip()
	return(R)
def Check_vcf2(Count):
	command = """head -{Count} Filter_GWAS/Merged.vcf | tail -1""".format(Count = Count + 39)
	O= Popen(command, shell=True,stdout=PIPE, stderr=PIPE)
	Out = O.communicate()[0].decode('UTF-8')
	R = Out.split()[7].split("R2=")[1].rstrip()
	return(R)


Output_name = "Output/{Metabo}.HRC.{Study}.EA.{Gender}.{DATE}.txt".format(Metabo=Metabo_dic[Metabo],Study=Study,Gender=Gender_dic[Gender],DATE=DT)
header = "SNPID,CHR,POSITION,EFFECT_ALLELE,OTHER_ALLELE,STRAND,BETA,SE,PVALUE,EAF,INFORMATIVE_AC,HWE_PVALUE,CALL_RATE,N_INFORMATIVE,IMPUTATION".replace(",","\t")
with open(Output_name,"w") as O:
	O.write(header+"\n")
Counter = 0
with gzip.open(Input_file,"r") as Test_results:
	for Loci in Test_results:
		Loci =  Loci.decode("utf-8") 
		if Loci[0] == "#": continue
		if "CHROM" in Loci: continue
		F = Loci.split()
		SNPID = F[0] + ":"+ F[1]
		CHR = F[0]
		POSITION = F[1]
		IMPUTATION = Check_vcf2(Counter)
		EFFECT_ALLELE = F[3]
		OTHER_ALLELE = F[2]
		STRAND = "NA"
		BETA = F[-2]
		SE = "NA"
		PVALUE = F[-1]  
		EAF = F[5]
		INFORMATIVE_AC = F[6]
		HWE_PVALUE = F[8]
		CALL_RATE = F[7]
		N_INFORMATIVE = F[4]	
		ROW = "\t".join([SNPID,CHR, POSITION,EFFECT_ALLELE, OTHER_ALLELE, STRAND, BETA, SE, PVALUE, EAF, INFORMATIVE_AC, HWE_PVALUE, CALL_RATE, N_INFORMATIVE, IMPUTATION])
		Counter += 1
		with open(Output_name,"a") as O:
			O.write(ROW + "\n")			
	
Expected = """
Column header 	Description 
SNPID 	SNP identification number
CHR	Chromosome number
POSITION 	Physical position for the reference sequence 
EFFECT_ALLELE 	Coded allele to which effect has been estimated (e.g. for A/G SNP in which AA=0, AG=1, GG=2, coded allele is G) 
OTHER_ALLELE	The other allele 
STRAND 	+ or -, representing either the positive/forward strand or the negative/reverse strand of the human genome reference sequence; to clarify which strand the coded_all and noncoded_all are on 
BETA 	Effect size for effect allele from genotype-phenotype association, at least 5 decimal places -- “NA” if not available
SE 	Standard error of the effect size, to at least 5 decimal places -- “NA” if not available
PVALUE 	P value for the association of variant and phenotype (uncorrected for genomic control) -- “NA” if not available 
EAF 	Frequency of the effect allele -- “NA” if not available
INFORMATIVE_AC	Allele count of the effect allele
HWE_PVALUE 	Hardy-Weinberg equilibrium p-value 
CALL_RATE 	Genotyping call rate 
N_INFORMATIVE	Total sample with phenotype and genotype for SNP 
IMPUTATION	A value (range 0-1) corresponding to the imputation quality measure (Rsq or info)
"""

