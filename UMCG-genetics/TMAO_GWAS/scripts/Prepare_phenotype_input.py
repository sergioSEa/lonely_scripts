#fid iid fatid matid sex y1 y2 y3 y4

covariates = []
phenotypes = []
with open("Transformed_covariates.tsv") as F:
	for line in F:
		if "age" in line: continue
		l = line.rstrip().split("\t")
		l[6] = str(int(l[6])+1)
		C = l[0:5]
		C.append(l[-1])
		C.append(l[6])
		phenotypes.append(C)
		covariates.append(l[5:])


Output_pheno = "Output/Phenotype.phe"
Output_cov = "Output/Covariates.phe"

Cov = [" ".join(["fid","iid","fatid","matid","sex","age","gender","bmi"])]
for line in covariates:
	New_line = [line[-1].strip('"'), line[-1].strip('"'), "0","0",line[1],line[0],line[1],line[2]]
	Cov.append(" ".join(New_line))
Pheno = [" ".join(["fid","iid","fatid","matid","sex","TMAO","Choline","Betaine","y-butyrobetaine","L-Carnitine"])]
for line in phenotypes:
	New_line = [line[-2].strip('"'), line[-2].strip('"'), "0","0",line[-1],line[0],line[1],line[2],line[3],line[4]]
	Pheno.append(" ".join(New_line))

with open(Output_pheno, "w") as O:
	O.write("\n".join(Pheno))
with open(Output_cov, "w") as O:
	O.write("\n".join(Cov))
	
