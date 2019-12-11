
Family_Latinname = "/Users/Sergio/Documents/Thesis/Rate_correlation/Species-Latin_name.txt"
Order_Family = "/Users/Sergio/Documents/Thesis/Rate_correlation/Family-Order.tsv"

dic_all = {}
with open(Family_Latinname,"r") as F:
	for line in F:
		family = line.rstrip().split()[0]
		name = line.rstrip().split()[1]
		if family not in dic_all:
			dic_all[family] = []
		dic_all[family].append(name)


dic_all2 = {}
dic_semi = {}
with open(Order_Family) as F2:
	for line in F2:
		l = line.rstrip().split()	
		Order = l[0]
		family = l[1]
		for name in dic_all[family]:
			dic_all2[name] = [family, Order]
			dic_semi[name.split("_")[0]] = ["NA", Order]


with open("Metadata_sp.csv") as M:
	for line in M:
		if "X1" in line: continue
		l = line.rstrip().replace('"','').split(",")
		
		if l[3] != "Birds" and l[4] != "Birds": continue
			
		if l[1] in dic_all2:
			print("\t".join([l[1], dic_all2[l[1]][0], dic_all2[l[1]][1]]))
		else:
			if l[1].split("_")[0] in dic_semi:
				m =  l[1].split("_")[0]
				print("\t".join([l[1], dic_semi[m][0], dic_semi[m][1]]))
			else:
				print("\t".join([l[1], "NA", "NA"]))
