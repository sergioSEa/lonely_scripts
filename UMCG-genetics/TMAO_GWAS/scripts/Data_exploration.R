library(tidyverse)


#Check distribution of Gender, BMI and Age
covariates = read_tsv("Data_Tables/covariates.tsv")
#age     gender_f1_m0    bmi


covariates %>% group_by(gender_f1_m0) %>% summary(n) -> SUMMARY_gender
covariates %>% drop_na() %>% select(-gender_f1_m0) %>% summarise_if(is.numeric, list(max, min, mean, median, sd)) -> SUMMARY
covariates %>% summarise_if(is.numeric, list(max, min, mean, median, sd)) -> SUMMARY


Make_manhattan = function(DF,NAME){
	GWAS = as_tibble(read.table(DF))
	colnames(GWAS) = as_vector(GWAS[1,])
	GWAS = GWAS[-1,]
	GWAS %>% filter(as.numeric(levels(IMPUTATION))[IMPUTATION] > 0.3) %>% filter(as.numeric(levels(EAF))[EAF] > 0.1) -> GWAS
	
	GWAS %>% mutate(PVALUE = as.numeric(levels(PVALUE))[PVALUE]) -> GWAS
	GWAS %>% arrange(as.numeric(levels(CHR))[CHR]) -> GWAS
	GWAS %>% mutate(Log_pvalue = -log10(as.numeric(PVALUE)), POSITION=as.numeric(POSITION), CS = cumsum(POSITION)) -> GWAS
	

	axis.set <- GWAS %>% group_by(CHR) %>% summarize(center = (max(CS) + min(CS)) / 2)
	ylim <- abs(floor(log10(min(as.numeric(GWAS$PVALUE))))) + 2   
	nCHR <- length(unique(GWAS$CHR))

	ggplot(GWAS, aes(x = CS, y = Log_pvalue, color = as.factor(CHR), size = Log_pvalue)) + geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(5e-8), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log10(p)") + theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )  -> GWAS_plot

ggsave(filename = paste(c("GWAS_plot", NAME, ".pdf"), collapse=""), GWAS_plot)



}


path <- "Output/"
file.names <- dir(path, pattern ="2020130.txt")
for(i in 1:length(file.names)){
	PATH = paste(c(path,file.names[i]), collapse="")
	print(PATH)
	Make_manhattan(PATH, file.names[i])
}






