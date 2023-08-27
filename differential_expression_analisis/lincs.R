Args <- commandArgs(T)
raw_expr <- Args[1]#the full path of gene expression profiles from LINCS (the output of get_LINCS_level3_gene_expr.py)
pheno_raw <- Args[2]#the full path of sample information (the output of get_LINCS_level3_gene_expr.py)
output <- Args[3]# the path of the output folder
setwd(output)

library(limma)
library(stringr)

data <- read.table(raw_expr,header=T,sep="\t",row.names = "pr_gene_symbol",check.names=FALSE)
pheno = read.table(pheno_raw,sep="\t",header=T)
pheno["drug_dose"] = paste(pheno$pert_iname,pheno$pert_dose,sep="_")
pheno["batch"] = paste(str_split_fixed(pheno$rna_plate,"_",n=4)[,1],str_split_fixed(pheno$rna_plate,"_",n=4)[,2],str_split_fixed(pheno$rna_plate,"_",n=4)[,3],sep="_")

all_drugs = unique(pheno$drug_dose)
print(paste("There are",length(all_drugs),"drug-dose pair",sep=" "))
print(all_drugs)


for (b in unique(pheno$batch)){
  print("-----------------------")
  print(b)
  dir.create(b)
  
  i=1
  for (drug in all_drugs[1:length(all_drugs)]){
    if (!str_detect(drug,"DMSO")){
      drug_expr = data[,as.character(pheno$inst_id[pheno$drug_dose==drug & pheno$batch == b])]
      
      if (class(drug_expr)!= "numeric"){
        if (dim(drug_expr)[2] >0 ) {
          print(drug)
          
          ctrl_id = pheno$inst_id[which(pheno$batch == b & pheno$pert_type=="ctl_vehicle")]
          ctrl_expr = data[,as.character(ctrl_id)]
          expr <- cbind(ctrl_expr,drug_expr)

          group <- c(rep("contrl",dim(ctrl_expr)[2]),rep("treatment",dim(drug_expr)[2]))
          f <- factor(group,levels = c("contrl","treatment"))
          
          design <- model.matrix(~0+f)
          colnames(design) <- levels(f)
          fit <- lmFit(expr,design)
          
          contrast <- makeContrasts(treatment-contrl,levels=design)
          fit <- contrasts.fit(fit,contrast)
          fit <- eBayes(fit)
          allDiff <- topTable(fit,adjust = 'BH',coef = 1,n = Inf)
          
          all_fc = as.data.frame(allDiff[,1],col.names=drug,row.names = row.names(allDiff))
          colnames(all_fc)[1] <- drug
          write.table(all_fc,file=paste(output,b,'/',drug,"_",b,".txt",sep=""),sep="\t",quote=F)
          
          print(i)
          i=i+1
        } 
        
      }
    } 
    
    
  }
}
