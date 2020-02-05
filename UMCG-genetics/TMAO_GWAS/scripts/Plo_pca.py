#import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

Data = pd.read_csv("plink.eigenvec", sep=" ")


print(Data.shape)


PC1 = Data.iloc[:,2]
PC2 = Data.iloc[:,3]


#ax = sns.scatterplot(x=PC1, y=PC2)
#ax.savefig("PCA.png")
plt.scatter(x=PC1, y=PC2)
plt.savefig(sys.argv[1])
plt.close


