from pathlib import Path
from Bio import SeqIO

def find_taxa(Accession):
	P = ""
	found = False
	taxonomy = "unknwon"
	for entry in Path("NCBI_SEQ/gbs").glob(Accession + ".[1-9].gb"):
		found = True
		with open(entry) as handle:
			S = SeqIO.parse(handle, "genbank")
			for entry in S:
				Tx = entry.annotations['taxonomy']
				if "Bifurcata" in Tx: taxonomy = "Bifurcata"
				elif "Aves" in Tx: taxonomy = "Aves"
				elif "Crocodylia" in Tx: taxonomy = "Crocodylia"
				elif "Testudines" in Tx: taxonomy = "Testudines"
				elif "Squamata" in Tx: taxonomy = "Squamata"
				elif "Mammalia" in Tx: taxonomy = "Mammalia"
				elif "Amphibia" in Tx: taxonomy = "Amphibia"
				elif "Actinopterygii" in Tx: taxonomy = "Actinopterygii"
				elif "Elasmobranchii" in Tx: taxonomy = "Elasmobranchii"
				elif "Chondrichthyes" in Tx: taxonomy = "Chondrichthyes"
				elif "Dipnoi" in Tx: taxonomy = "Dipnoi"
				else:
					try: 
						p = Tx.index("Vertebrata")
						taxonomy = Tx[p+1]
					except:				
						taxonomy = "unknown"

	if found == True:
		return(taxonomy)	

	for entry in Path("NCBI_SEQ/complete_gbs").glob(Accession + ".*.gb"):
		with open(entry) as handle:
			S = SeqIO.parse(handle, "genbank")
			for entry in S:
				Tx = entry.annotations['taxonomy']
				if "Bifurcata" in Tx: taxonomy = "Bifurcata"
				elif "Aves" in Tx: taxonomy = "Aves"
				elif "Crocodylia" in Tx: taxonomy = "Crocodylia"
				elif "Testudines" in Tx: taxonomy = "Testudines"
				elif "Squamata" in Tx: taxonomy = "Squamata"
				elif "Mammalia" in Tx: taxonomy = "Mammalia"
				elif "Amphibia" in Tx: taxonomy = "Amphibia"
				elif "Actinopterygii" in Tx: taxonomy = "Actinopterygii"
				elif "Elasmobranchii" in Tx: taxonomy = "Elasmobranchii"
				elif "Chondrichthyes" in Tx: taxonomy = "Chondrichthyes"
				elif "Dipnoi" in Tx: taxonomy = "Dipnoi"
				else: 
					try: 
						p = Tx.index("Vertebrata")
						taxonomy = Tx[p+1]
					except:
						taxonomy = "unknown"
	return(taxonomy)



with open("ND3_genes.fa") as F:
#with open("ND3_complete_genes.fa") as F:
	for line in F:
		if line[0] ==">":
			acc = line.rstrip().split("-")[1]	
			tax = find_taxa(acc)
			print(line[1:-1],tax)
