
#devtools::install_github("YuLab-SMU/ggtree")
#devtools::install_github("ssun6/plotmicrobiome")

library(plotmicrobiome)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)

setwd("/Users/shansun/git/MicroscopicColitisMicrobiome")
metadata_dir="./metadata/metadata_short.csv"
taxa_dir="./data/taxa_combined.csv"

#format ASV table and metadata
taxa_tab1=format_asv(taxa_file = taxa_dir,biom=F,onefile = T,ASV=F,sep=",")
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=1)

#Fig.1
#PCoA plot and alpha diversity
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="location",stratify_by_value="ASC",prevalence_cutoff=0, abundance_cutoff=0)
pdf("./output/pcoa_asc.pdf",height=15,width=5,onefile=T)
mds_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="Case_Ctrl",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")
dev.off()
pdf("./output/alpha_asc.pdf",height=30,width=10,onefile=T)
alpha_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="Case_Ctrl",palette_group=c("red","blue","orange","green"))
dev.off()

tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="location",stratify_by_value="DES",prevalence_cutoff=0, abundance_cutoff=0)
pdf("./output/pcoa_des.pdf",height=15,width=5,onefile=T)
mds_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="Case_Ctrl",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")
dev.off()
pdf("./output/alpha_des.pdf",height=30,width=10,onefile=T)
alpha_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="Case_Ctrl",palette_group=c("red","blue","orange","green"))
dev.off()

#model results
#ASC
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="location",stratify_by_value="ASC",prevalence_cutoff=0, abundance_cutoff=0)
tax_l=sapply(strsplit(rownames(tab_s),"--"),length)
tab_s1=tab_s[which(tax_l==6),]#genus level
gen_mds=vegan::capscale(t(tab_s1)~1,distance="bray")
mat1=summary(gen_mds)$sites
alpha1=diversity(tab_s1,index="shannon",MARGIN = 2)
mat2=cbind(alpha1,mat1)
metadata2=metadata1[match(rownames(mat2),rownames(metadata1)),]
colnames(mat2)=c("alpha_ASC","MDS1_ASC","MDS2_ASC","MDS3_ASC","MDS4_ASC","MDS5_ASC","MDS6_ASC")

models=list()
models[[1]]="Case_Ctrl+factor(grade_school)+factor(ppi)+factor(batch)"
models[[2]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic"
models[[3]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic+patient_age"
models[[4]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic+patient_age+BMI"

bacteriaMeta1=cbind(mat2,metadata2)
pVals=matrix(ncol=7,nrow=4)
tVals=matrix(ncol=7,nrow=4)
for (m in 1:4){
  for (i in 1:7){
    model <- as.formula(paste(colnames(bacteriaMeta1)[i],"~",models[[m]]))
    
    simpleMod <- try(lm(model,data=bacteriaMeta1))
    
    if(class( simpleMod )=="try-error"){
      pVals[i,]  <-NA
      tVals[i,]  <-NA
    }
    pVals[m,i] <- summary(simpleMod)$coefficients[2,4]
    tVals[m,i] <- summary(simpleMod)$coefficients[2,3]
  }
}
mds_vals=cbind(tVals,pVals)
rownames(mds_vals)=c("model1","model2","model3","model4")
colnames(mds_vals)=paste(c("alpha_ASC","MDS1_ASC","MDS2_ASC","MDS3_ASC","MDS4_ASC","MDS5_ASC","MDS6_ASC"),
                         c(rep("t",7),rep("P",7)),sep="_")
write.csv(mds_vals,file="./output/mds_models_asc.csv")

#DES
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="location",stratify_by_value="DES",prevalence_cutoff=0, abundance_cutoff=0)
tax_l=sapply(strsplit(rownames(tab_s),"--"),length)
tab_s1=tab_s[which(tax_l==6),]#genus level
gen_mds=vegan::capscale(t(tab_s1)~1,distance="bray")
mat1=summary(gen_mds)$sites
alpha1=diversity(tab_s1,index="shannon",MARGIN = 2)
mat2=cbind(alpha1,mat1)
metadata2=metadata1[match(rownames(mat2),rownames(metadata1)),]
colnames(mat2)=c("alpha_DES","MDS1_DES","MDS2_DES","MDS3_DES","MDS4_DES","MDS5_DES","MDS6_DES")

models=list()
models[[1]]="Case_Ctrl+factor(grade_school)+factor(ppi)+factor(batch)"
models[[2]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic"
models[[3]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic+BMI"
models[[4]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic+BMI+patient_age"

bacteriaMeta1=cbind(mat2,metadata2)
pVals=matrix(ncol=7,nrow=4)
tVals=matrix(ncol=7,nrow=4)
for (m in 1:4){
  for (i in 1:7){
    model <- as.formula(paste(colnames(bacteriaMeta1)[i],"~",models[[m]]))
    
    simpleMod <- try(lm(model,data=bacteriaMeta1))
    
    if(class( simpleMod )=="try-error"){
      pVals[i,]  <-NA
      tVals[i,]  <-NA
    }
    pVals[m,i] <- summary(simpleMod)$coefficients[2,4]
    tVals[m,i] <- summary(simpleMod)$coefficients[2,3]
  }
}
mds_vals=cbind(tVals,pVals)
rownames(mds_vals)=c("model1","model2","model3","model4")
colnames(mds_vals)=paste(c("alpha_DES","MDS1_DES","MDS2_DES","MDS3_DES","MDS4_DES","MDS5_DES","MDS6_DES"),
                         c(rep("t",7),rep("P",7)),sep="_")
write.csv(mds_vals,file="./output/mds_models_des.csv")


#Fig. 2
#subset the abundance table to only include samples for test
taxa_tab1=format_asv(taxa_file = taxa_dir,biom=F,onefile = T,ASV=F,sep=",")
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=1)
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="",stratify_by_value="",prevalence_cutoff=0.1, abundance_cutoff=0)
tax_l=sapply(strsplit(rownames(tab_s),"--"),length)
tab_s=tab_s[which(tax_l!=8),]#genus level
tab1=tab_s[,match(intersect(colnames(tab_s),rownames(metadata1)),colnames(tab_s))]
metadata1=metadata1[match(intersect(colnames(tab_s),rownames(metadata1)),rownames(metadata1)),]

models=list()
models[[1]]="Case_Ctrl+factor(grade_school)+factor(ppi)+factor(batch)"
models[[2]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic"
models[[3]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic+patient_age"
models[[4]]="Case_Ctrl+Sex+factor(grade_school)+factor(ppi)+factor(batch)+antibiotic+patient_age+BMI"

taxa_table5=tab1
loc1=c("ASC","DES")
rownames(taxa_table5)=paste0("a",c(1:nrow(taxa_table5)))
bacteriaMeta1=cbind(t(taxa_table5),metadata1)
pVals=matrix(ncol=8,nrow=nrow(taxa_table5))
tVals=matrix(ncol=8,nrow=nrow(taxa_table5))
for (m in 1:4){
  for (j in 1:2){
    bacteriaMeta=bacteriaMeta1[bacteriaMeta1$location==loc1[j],]
    for( i in 1:dim(taxa_table5)[1]){
      print(i)
      if(sd(taxa_table5[i,])==0){
        next
      }
      model <- as.formula(paste(rownames(taxa_table5)[i],"~",models[[m]]))
      
      simpleMod <- try(lm(model,data=bacteriaMeta,na.action=na.omit))
      
      if(class( simpleMod )=="try-error"){
        pVals[i,]  <-NA
        tVals[i,]  <-NA
      }
      pVals[i,(m-1)*2+j] <- summary(simpleMod)$coefficients[2,4]
      tVals[i,(m-1)*2+j] <- summary(simpleMod)$coefficients[2,3]
    }
  }
}

fdrs=pVals
for(n in 1:8){
  fdrs[,n]=p.adjust(pVals[,n],method="fdr")
}

rownames(tVals)=rownames(tab_s)
rownames(pVals)=rownames(tab_s)
rownames(fdrs)=rownames(tab_s)

colnames(tVals)=c("model1_ASC","model1_DES","model2_ASC","model2_DES","model3_ASC","model3_DES","model4_ASC","model4_DES")
colnames(pVals)=c("model1_ASC","model1_DES","model2_ASC","model2_DES","model3_ASC","model3_DES","model4_ASC","model4_DES")
colnames(fdrs)=c("model1_ASC","model1_DES","model2_ASC","model2_DES","model3_ASC","model3_DES","model4_ASC","model4_DES")

fdrs_min=apply(fdrs,1,min)
which(fdrs<0.1,T)
colnames(pVals)[which(fdrs<0.1,T)[,2]]

length(which(fdrs[,1]<0.1))
length(which(fdrs[,2]<0.1))
length(which(fdrs[,3]<0.1))
length(which(fdrs[,4]<0.1))
length(which(fdrs[,5]<0.1))
length(which(fdrs[,6]<0.1))
length(which(fdrs[,7]<0.1))
length(which(fdrs[,8]<0.1))

length(which(fdrs[,2]<0.1 & fdrs[,4]<0.1))
length(which(fdrs[,2]<0.1 & !fdrs[,4]<0.1))
length(which(!fdrs[,2]<0.1 & fdrs[,4]<0.1))

length(which(fdrs[,2]<0.1 & tVals[,2]<0))
length(which(fdrs[,2]<0.1 & tVals[,2]>0))

length(which(fdrs[,4]<0.1 & tVals[,4]<0))
length(which(fdrs[,4]<0.1 & tVals[,4]>0))


write.csv(pVals,file="./output/lmer_P_cov_all_model.csv")
write.csv(fdrs,file="./output/lmer_FDR_cov_all_model.csv")

r_all=cbind(tVals,pVals,fdrs)
colnames(r_all)=paste(colnames(r_all),c(rep("t",8),rep("P",8),rep("FDR",8)),sep="_")
r_all=r_all[order(r_all[,2]),]
write.csv(r_all,file="./output/lmer_all_models_results.csv")

r_all_sig=r_all[which(apply(r_all[,17:24],1,min)<0.1),]
write.csv(r_all_sig,file="./output/lmer_all_models_sig_results.csv")


#Fig.2
pdf("./output/tree_model.pdf",height =5, width=10,onefile=T)
for (i in 1:4){
  #only DES trees because there is no significant ones in ASC
  fdrs1=cbind(tVals[,i*2],pVals[,i*2],fdrs[,i*2])
  rownames(fdrs1)=rownames(fdrs)
  head(fdrs1[order(fdrs1[,1]),],n=10)
  
  tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="location",stratify_by_value="DES",prevalence_cutoff=0, abundance_cutoff=0)
  tab1=tab_s[,intersect(colnames(tab_s),rownames(metadata1))]
  metadata1=metadata1[match(intersect(colnames(tab_s),rownames(metadata1)),rownames(metadata1)),]
  
  plot1=tree_view(taxa_table =tab1, metadata=metadata1,fdrs=fdrs1,test_metadata="Case_Ctrl",prevalence_cutoff=0.05, abundance_cutoff=0)
  print(plot1)
}
dev.off()


#Fig.3
pdf("./output/corplot_model.pdf",height =8, width=8,onefile=T)
pVals_log=data.frame(-log10(fdrs)*sign(tVals))
cor1=matrix(nrow=4,ncol=2)
for (i in 1:4){
  p=ggplot(pVals_log, mapping=aes_string(x=colnames(pVals_log)[i*2-1], y=colnames(pVals_log)[i*2])) +geom_point(color = 'red')+ theme_classic(base_size = 20) + labs(title="Comparison of ASC and DES",x ="ASC -log10(FDR)*direction" , y = "DES -log10(FDR)*direction")
  tax_lab=rownames(pVals_log)
  tax_lab1=sapply(strsplit(tax_lab,"--"),function(i){i[length(i)]})
  tax_lab2=sapply(strsplit(tax_lab,"--"),function(i){i[length(i)-1]})
  tax_lab1[which(tax_lab1=="__")]=paste0(tax_lab2[which(tax_lab1=="__")],"__unclassified")
  
  tax_lab1[which(fdrs[,i*2-1]>0.1 & fdrs[,i*2]>0.1)]=NA
  p2=p+geom_text_repel(aes(label =tax_lab1),size = 4)+geom_vline(xintercept=0, linetype="dotted")+geom_hline(yintercept=0, linetype="dotted")
  print(p2)
  cor1[i,1]=cor.test(pVals_log[,i*2-1],pVals_log[,i*2],method="spearman")$estimate
  cor1[i,2]=cor.test(pVals_log[,i*2-1],pVals_log[,i*2],method="spearman")$p.value
}
dev.off()
colnames(cor1)=c("rho","P")
rownames(cor1)=c("model1","model2","model3","model4")
write.csv(cor1,file="./output/corplot_model.csv")
