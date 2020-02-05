

read.table("Data_Tables/Total_merged_Covariates.tsv")


Covariates = read.table("Data_Tables/Total_merged_Covariates.tsv", sep="\t", header=T)

Covariates = Covariates[!is.na(Covariates$Age..years.),]
Covariates = Covariates[!is.na(Covariates$Sex),]
Covariates = Covariates[!is.na(Covariates$bmi),]
Covariates[complete.cases(Covariates[ , c("bmi","Sex","Age..years.","TMAO")]),] -> Covariates

Model = lm(TMAO ~ `Age..years.` + Sex + bmi, Covariates)


DATA = data.frame(Covariates$LLDEEPID.1, Model$residuals)
x = DATA[,2]
y<-qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))

DATA[,2] = y

write.table(DATA, file='covariate_corrected.tsv', quote=FALSE, sep='\t', col.names = NA)

