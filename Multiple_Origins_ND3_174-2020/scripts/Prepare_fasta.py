#Get fasta file for all gbf files
#Create taxonomy file in the process


from pathlib import Path
from Bio import SeqIO
#Create a new Taxonomy file with information about the taxonomy for each gbf records
with open("taxonomy", "w") as O:
	pass	

MISS = []
Output = ""


for File in Path("NCBI_SEQ/complete_gbs").glob("*.gb"):
#for File in Path("../Gb2Tree_reptile_and_bird_mtgenomes/").glob("*/gbs/*.gb"):
	IN = False
	with open(File,"r") as handle:
		Entry = SeqIO.parse(handle,"genbank")
		for i in Entry:
			NCBI_ID = i.name
			organism = "_".join(i.description.split()[0:2])
			sequence = str(i.seq)
			#Taxonomy assignment of the sequence
			Tx = i.annotations['taxonomy']
			if "Bifurcata" in Tx: taxonomy = "Bifurcata"
			elif "Aves" in Tx: taxonomy = "Aves"
			elif "Crocodylia" in Tx: taxonomy = "Crocodylia"
			elif "Testudines" in Tx: taxonomy = "Testudines"
			elif "Squamata" in Tx: taxonomy = "Squamata"
			elif "Mammalia" in Tx: taxonomy = "Mammalia"
			elif "Amphibia" in Tx: taxonomy = "Amphibia"
			elif "Actinopterygii" in Tx: taxonomy = "Actinopterygii"
			elif "Elasmobranchii" in Tx: taxonomy = "Elasmobranchii"
			else: taxonomy = "unknown" #Records with unkown category will be checked by hand.
			#Look for ND3 record. Gbf files might be different.
			for gene in i.features:
				try: g = gene.qualifiers["gene"][0].strip("'")
				except: 
					try:
						if gene.type == "CDS":
							g = gene.qualifiers["product"][0].strip('"')
						else: continue
					except:
						continue
				if g == "ND3" or  g== "NADH dehydrogenase subunit 3":
					IN = True
					sequence = str(i.seq)			
					l = list(gene.location)
					try: b= l[0]
					except: continue
					e=l[-1]
					s = sequence[b:e]
		#If ND3 was found in the record, save taxonomy and print the record out
		if IN == True:		
			Output = ">" + organism + "-" + NCBI_ID + "\n" + s
			with open("taxonomy", "a") as O:
				O.write(organism + "\t" + taxonomy + "\n")
			print(Output)
		else:
			MISS.append(str(File))


#Save the records that are missing the ND3 sequence
with open("NCBI_SEQ/gbs/Missing_ND3.txt","w") as F:
	F.write("\n".join(MISS))			
