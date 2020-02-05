to_remove = []
with open("Output/plink.genome") as Info_table:
	for line in Info_table:
		try : Pi= float(line.split()[9])
		except: continue
		if Pi >= 0.2: 
			to_remove.append(line.rstrip())

with open("Output/related.tsv","w") as F:
	F.write("\n".join(to_remove))
for i in to_remove:
	
	ID = "_".join(i.split()[0:2])
	ID2 = "_".join(i.split()[2:4])
	print(ID)
	#print(ID2)
