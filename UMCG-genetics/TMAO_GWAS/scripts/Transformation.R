#TMAO    Choline Betaine y-butyrobetaine L-Carnitine     age     gender  bmi     ID
args = commandArgs(trailingOnly=TRUE)

Phenotypes = read.table(args[1], sep="\t", header=T)

Phenotypes[complete.cases(Phenotypes[ , c("bmi","gender","age")]),] -> Phenotypes


Transform = function(V){
y<-qnorm((rank(V,na.last="keep")-0.5)/sum(!is.na(V)))
return(y)
}


for (i in c("TMAO","Choline","Betaine","y.butyrobetaine","L.Carnitine")){
	Number = which(colnames(Phenotypes) == i)
	Phenotypes[,Number] = Transform(Phenotypes[,Number])
}

write.table(Phenotypes,"Transformed_covariates.tsv",sep="\t",row.names=F)
