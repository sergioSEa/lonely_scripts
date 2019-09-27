import numpy as np
from pathlib import Path
import yaml
from subprocess import call

with open("config.yaml", 'r') as stream:
	try:
		input_dic = yaml.safe_load(stream)
	except yaml.YAMLError as exc:
		print(exc)


FILE = input_dic['SNV_call']
path_coverage1 =  input_dic['average_coverage']
path_coverage2 = input_dic['percentage_coverage']
relatedness = input_dic['relation'] #Vector specifying which columns of the variant call should be added together
outdir= input_dic["output_dir"]
VC = input_dic["mean_depth_filter"]
HC =  input_dic["coverage_filter"]

dic_vertical = {}
dic_horizontal={} 

with open(path_coverage1) as File:
	for line in File:
		if "Sample" in line: continue
		l = line.split()
		sc = l[0]
		vertical = l[1:]
		dic_vertical[sc] = vertical
with open(path_coverage2) as File:
        for line in File:
                if ".bam" in line: continue
                l = line.split()
                sc = l[0]
                horizontal = l[1:]
                dic_horizontal[sc] = horizontal
##Check if number of samples given in relatedness and the number of samples in the files match
for item in dic_horizontal.keys():
	if len(relatedness) != len(dic_horizontal[item]):
		expected =  len(dic_horizontal[item])
		given = len(relatedness)
		exit("Error: Number of samples given in 'Relation' does not match the number of samples being analyzed. Expected: {Ex}, Given: {Giv}".format(Ex=expected, Giv=given))
	else: break

if not Path(outdir).exists():
	print("Creating output directory:" + outdir )
	command = "mkdir "+outdir
	call(command,shell=True)


missing= outdir + "/filter_SNP.txt"
output = outdir + "/frequency_SNP.txt"



with open(missing,"w") as F: 
	pass

def create_vector(Input,depth,coverage,vec=relatedness):
	#FILTER BY VERTICAL COVERAGE (MEAN DEPTH OF A SCAFFOLD) --> VC
	#FILTER BY HORIZONTAL COVERAGE (COVERAGE OF THE SCAFFOLD) --> HC
	splitted= Input.split("|")

	for n in range(len(coverage)):
		#The length of coverage is the number of samples. Iteration to the coverage and mean depth of coverage of each of the samples and check if the sample should be filtered out.
		c = coverage[n]
		d = depth[n]
		if c < HC or d < VC:
			splitted[n] = 0
	splitted= [ float(x) for x in splitted ]
	#Taking into account vector of relations:
	dic_elements = {}
	counter = 0 
	for element in vec:
		if element not in dic_elements:
			dic_elements[element] = []
		dic_elements[element].append(counter)
		counter += 1
	splitts = []
	for element in dic_elements:
		list_positions = dic_elements[element]
		total = 0
		for i in list_positions:
			total += splitted[i]
		splitts.append(total)
			
	spllits= np.array(splitts)

	return(spllits)


def Iteration():
	MATRIX_snp = None
	np.seterr(divide='ignore', invalid='ignore')
	print("Starting filtering of Variant calls")
	with open(FILE) as F:
		#Big variant call file with a SNP per line, info on which scaffold is present, total counts and counts that are a non-reference allele
		for line in F:
			l = line.split()
			scaffold = l[0]
			pos= l[2]
			ref= l[3] #reference allele
			#### Line-separed fields
			total_counts=l[4]
			alt_counts=l[5]
			####
			name= scaffold +"_"+ pos + "_"+ ref #unique identifier
	
			# Filters
			horizontal_filter = dic_horizontal[scaffold] #coverage of the scaffold
			vertical_filter = dic_vertical[scaffold] #mean depth of coverage of the scaffold
		
			horizontal_filter = [ float(x) for x in horizontal_filter ]
			vertical_filter = [ float(x) for x in vertical_filter ]
		
			# Uses the relation vector given to generate the total number of counts, filters according to the filters given		
			counts = create_vector(total_counts,vertical_filter,horizontal_filter)
			#Get samples where there's no counts in that position
			missi = counts[counts==0]
		
			if len(missi) >= len(relatedness) -1:
				with open(missing, "a") as R:
					R.write(line)
					continue
			if "," in alt_counts:
				for variant in alt_counts.split(","):
					v = variant.split(".|")
					n_v = v[0].split("|")[1]
					variant_counts = create_vector(v[1],vertical_filter,horizontal_filter)
					np.seterr(divide='ignore')
					frequency = np.divide(variant_counts,counts)
			
					if np.all(frequency == 0): continue
					name +=  "_"+n_v
					
					try:
						if MATRIX_snp == None: MATRIX_snp = frequency
						else: MATRIX_snp = np.concatenate(MATRIX_snp,frequency)
					except: MATRIX_snp = np.vstack((MATRIX_snp,frequency))
					with open(output,"w") as G:
						G.write(name +" "+ str(frequency[0]) +" "+ str(frequency[1]))

			else:
				v = alt_counts.split(".|")
				n_v = v[0].split("|")[1]
				variant_counts = create_vector(v[1],vertical_filter,horizontal_filter)
				frequency = variant_counts/counts
				if np.all(frequency == 0): continue
				name +=  "_"+n_v
				#print(name,frequency)
				try:
					if MATRIX_snp == None: MATRIX_snp = frequency
					else: MATRIX_snp = np.concatenate(MATRIX_snp,frequency)
				except:  MATRIX_snp = np.vstack((MATRIX_snp,frequency))
				with open(output,"w") as G:
					 G.write(name +" "+ str(frequency[0]) +" "+ str(frequency[1]))

	return(MATRIX_snp)
def distance(x,y):
	return(sum((x - y)**2)**0.5)

def pi(v1,v2):
	F1 = 1-v1
	F2 = 1-v2
	
	PI = v1*F2 + v2*F1
	PI = sum(PI) / len(v1)
	return(PI)
	
def Fst(v1,v2):
	if len(v1) == 0: 
		Fst = np.nan
		return(Fst)
	numerator = (pi(v1,v1) + pi(v2,v2))/2
	denominator = pi(v1,v2)
	Fst = 1 - (numerator/denominator)
	return(Fst)

from scipy import spatial,stats

MATRIX_snp = Iteration()
X = MATRIX_snp.transpose()

#X = np.array([[0.2,0.1,0.7], [0.2,0.5,np.nan], [0,0,0.8]])



#np.save("Matrix",X)


#Iterete through organism and calculate distance between them in the positions both have information
#Row = organism ; column = SNP

distance = np.zeros((X.shape[0],X.shape[0]))
Fst_matrix =  np.zeros((X.shape[0],X.shape[0]))
print("Distance computation")
for item in range(X.shape[0]):
	for  item2 in range(X.shape[0]):
		x = X[item]
		y = X[item2]
		x = x.reshape(1,x.shape[0])
		y = y.reshape(1,y.shape[0])
	
		Y = np.concatenate((x,y))
		Y = Y.transpose()	
		Y = Y[~np.isnan(Y).any(axis=1)]
		Y = Y[np.isfinite(Y).any(axis=1)]
		Y = Y.transpose()
		if len(Y[0]) == 0:
			D = np.nan
		else:
			D = spatial.distance.pdist(Y, metric='cityblock')
			if len(D) == 0: D = np.nan
			else: D = D[0]
		
		FST = Fst(Y[0],Y[1])

		Fst_matrix[item,item2] = FST
		distance[item, item2] = D
print(Fst_matrix)
print(distance)		


output2 = outdir + "/manhattan_dist.txt"
output3 =  outdir + "/Fst_dist.txt"
np.savetxt(output2, distance, delimiter="\t")
np.savetxt(output3, Fst_matrix, delimiter="\t")



exit()


#X = stats.zscore(X)

#print(X)

#exit()


eu = spatial.distance.pdist(X, metric='euclidean')
manh  =  spatial.distance.pdist(X, metric='cityblock')

print(distance(X[0,],X[1,]))
print(eu)
print(manh)

