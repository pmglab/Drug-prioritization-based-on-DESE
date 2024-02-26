

Args <- commandArgs(T)
raw_expr <- Args[1] #drug and vehicle-induced gene expression profiles generated from get_LINCS_level3_gene_expr.py
pheno_raw <- Args[2]#drug and vehicle-induced gene expression profiles-related information generated from get_LINCS_level3_gene_expr.py
output <- Args[3]# folder to place the differential expression analysis

library(limma)
library(stringr)

setwd(output)

#load gene expression data
data <- read.table(raw_expr,header=T,sep="\t",row.names = "pr_gene_symbol",check.names=FALSE)
#load drug information
pheno = read.table(pheno_raw,sep="\t",header=T,quote="")

#'drug-dose' pair
pheno["drug_dose"] = paste(pheno$pert_iname,pheno$pert_dose,sep="_")

#plate type
pheno["plate"] = paste(str_split_fixed(pheno$rna_plate,"_",n=4)[,1],str_split_fixed(pheno$rna_plate,"_",n=4)[,2],str_split_fixed(pheno$rna_plate,"_",n=4)[,3],sep="_")

print(head(pheno))
#statistic the total 'drug-dose' pair
all_drugs = unique(pheno$drug_dose)
print(paste("There are",length(all_drugs),"drug-dose pair.",sep=" "))
print(all_drugs)

#compute the standard deviation of the log2foldchanges of a certain gene induced by a treatment
std.error <- function(x) sd(x)/sqrt(length(x))

#perform differential gene expression analysis based on the drug-induced gene expession profiles and vehicle-induced gene expression profiles on the same plate type
for (plt in unique(pheno$plate)){
  print("-----------------------")
  print(plt)
  dir.create(plt)
  
  i=1
  for (drug in all_drugs[1:length(all_drugs)]){
    if (!str_detect(drug,"DMSO")){
      
      #extract drug-induced gene expression profiles
      drug_expr = data[,as.character(pheno$inst_id[pheno$drug_dose==drug & pheno$plate == plt])]
      if (class(drug_expr)!= "numeric"){
        if (dim(drug_expr)[2] >0 ) {
          print(drug)
          
          #extract control vehicle-induced gene expression profiles
          ctrl_id = pheno$inst_id[which(pheno$plate == plt & pheno$pert_type=="ctl_vehicle")]
         
          if (length(ctrl_id)[1]>0){

          ctrl_expr = data[,as.character(ctrl_id)]

          #compute the se
          ctrl_expr_mean = apply(ctrl_expr,MARGIN=1,FUN=mean)
          fc = drug_expr - ctrl_expr_mean
          fc_se = apply(fc,MARGIN=1,FUN=std.error)

          #combine drug-induced and vehicle-induced gene expression profiles
          expr <- cbind(ctrl_expr,drug_expr)

          #differential expression analysis
          group <- c(rep("contrl",dim(ctrl_expr)[2]),rep("treatment",dim(drug_expr)[2]))
          f <- factor(group,levels = c("contrl","treatment"))
          design <- model.matrix(~0+f)
          colnames(design) <- levels(f)
          fit <- lmFit(expr,design)
          contrast <- makeContrasts(treatment-contrl,levels=design)
          fit <- contrasts.fit(fit,contrast)
          fit <- eBayes(fit)
          allDiff <- topTable(fit,adjust = 'BH',coef = 1,n = Inf)
          
          #keep the log2FoldChanges values and se values
          all_fc = as.data.frame(allDiff[,1],col.names=drug,row.names = row.names(allDiff))
          colnames(all_fc)[1] <- drug
          all_fc = merge(all_fc,fc_se,by.x=0,by.y=0)
          colnames(all_fc)[3] <- paste(drug,".SE",sep="")
          write.table(all_fc,file=paste(output,plt,'/',drug,"_",plt,".txt",sep=""),sep="\t",quote=F,row.names=F)
          
          print(i)
          i=i+1
        }
        } 
        
      }
    } 
    
    
  }
}
