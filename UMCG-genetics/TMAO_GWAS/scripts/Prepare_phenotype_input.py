#fid iid fatid matid sex y1 y2 y3 y4

covariates = []
phenotypes = []
with open("Transformed_covariates.tsv") as F:
	for line in F:
		#['"ID"', '"TMAO"', '"Choline"', '"Betaine"', '"y.butyrobetaine"', '"L.Carnitine"', '"TMAO.Choline"', '"TMAO.Betaine"', '"TMAO.Butyrobetaine"', '"TMAO.Carnitine"', '"Butyrobetain.Carnitine"', '"Age"', '"Gender"', '"BMI"']
		if "Age" in line: continue
		l = line.rstrip().split("\t")
		l[0] = l[0].strip('"') #"1_"+l[0].strip('"')
		l[-2] = str(int(l[-2])+1)
		C = l[1:11]
		C.append(l[0])
		C.append(l[-2])
		phenotypes.append(C)
		covariates.append( [l[0],l[-2], l[-3], l[-1] ]  )
		

Output_pheno = "Output/Phenotype.phe"
Output_cov = "Output/Covariates.phe"

Cov = [" ".join(["fid","iid","fatid","matid","sex","age","gender","bmi"])]
for line in covariates:
	New_line = [line[0].strip('"'), line[0].strip('"'), "0","0",line[1],line[2],line[1],line[3]]
	Cov.append(" ".join(New_line))
Pheno = [" ".join(["fid","iid","fatid","matid","sex","TMAO","Choline","Betaine","y-butyrobetaine","L-Carnitine","TMAO.Choline", "TMAO.Betaine", "TMAO.Butyrobetaine", "TMAO.Carnitine", "Butyrobetaine.Carnitine"])]
for line in phenotypes:
	New_line = [line[-2].strip('"'), line[-2].strip('"'), "0","0",line[-1],line[0],line[1],line[2],line[3],line[4], line[5], line[6], line[7], line[8], line[9]]
	Pheno.append(" ".join(New_line))
	
with open(Output_pheno, "w") as O:
	O.write("\n".join(Pheno))
with open(Output_cov, "w") as O:
	O.write("\n".join(Cov))
	
