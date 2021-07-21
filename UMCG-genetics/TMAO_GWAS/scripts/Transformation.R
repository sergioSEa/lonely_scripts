#TMAO    Choline Betaine y-butyrobetaine L-Carnitine     age     gender  bmi     ID
args = commandArgs(trailingOnly=TRUE)

Phenotypes = read.table(args[1], sep="\t", header=T)

Phenotypes[complete.cases(Phenotypes[ , c("BMI","Gender","Age")]),] -> Phenotypes


Transform = function(V){
y<-qnorm((rank(V,na.last="keep")-0.5)/sum(!is.na(V)))
return(y)
}


for (i in c("TMAO","Choline","Betaine","y.butyrobetaine","L.Carnitine", "TMAO.Choline", "TMAO.Betaine", "TMAO.Butyrobetaine", "TMAO.Carnitine", "Butyrobetain.Carnitine")){
	Number = which(colnames(Phenotypes) == i)
	Phenotypes[,Number] = Transform(Phenotypes[,Number])
}

write.table(Phenotypes,"Transformed_covariates.tsv",sep="\t",row.names=F)
