Args <- commandArgs(T)
raw_expr <- Args[1]#the output of get_GSE119291_expr.R; the full path of "all_hips_ctl_expression.txt" or "all_hips_sz_expression.txt" (the output of get_GSE119291_expr.R)
pheno_raw <- Args[2]#the full path of sample information getting from GEO website(GSE119291)
type <- Args[3]#str type, "hiPSC_control" or "hiPSC_SZ"
output <- Args[4]# the path of the output folder

library(limma)
library(stringr)

data = read.table(raw_expr,header=T,sep="\t",row.names = "Group.1",check.names=FALSE)
pheno = read.table(pheno_raw,sep="\t",header=T)
pheno = pheno[which(pheno$Source.name==type),]
pheno["drug_dose"] = paste(pheno$Perturbagen,pheno$Dosage,sep="_")
print(head(pheno))

all_drugs = unique(pheno$drug_dose)
print(paste("There are",length(all_drugs),"drug-dose pair",sep=" "))
print(all_drugs)

ctrl = data[,as.character(pheno[which(pheno$Perturbation.type=="vehicle"),]$Accession)]
cmm = merge(pheno,as.data.frame(t(ctrl)),by.x="Accession",by.y="row.names")
cmm1 = cmm[,c("Cell.id",colnames(t(ctrl)))]
cmm2 = aggregate(cmm1[,2:dim(cmm1)[2]],by=list(cmm1$Cell.id),FUN = mean)
cmm3 = as.data.frame(cmm2[,2:dim(cmm2)[2]],row.names = as.character(cmm2$Group.1))
ctrl_merged = t(cmm3)

i=1
for (drug in all_drugs[1:length(all_drugs)]){
  if (!str_detect(drug,"DMSO")){
    print(drug)
    drug_expr = data[,as.character(pheno$Accession[pheno$drug_dose==drug])]
    
    if (class(drug_expr)!= "numeric"){
        com_indiv = intersect(colnames(ctrl_merged),unique(as.character(pheno$Cell.id[pheno$drug_dose==drug])))
        print(com_indiv)
        
        print(paste(drug,"has >=2 batch.",sep=" ")) 
        drug1 = merge(pheno,t(drug_expr),by.x="Accession",by.y=0)
        drug2 = drug1[,c("Cell.id",colnames(t(drug_expr)))]
        drug3 = aggregate(drug2[,2:dim(drug2)[2]],by=list(drug2$Cell.id),FUN = mean)
        drug4 = as.data.frame(drug3[,2:dim(drug3)[2]],row.names = as.character(drug3$Group.1))
        drug_expr1 = t(drug4)
        
        ctrl_1 = ctrl_merged[,com_indiv]
        drug_expr1 = drug_expr1[,com_indiv]
        
        colnames(ctrl_1) <- paste(colnames(ctrl_1),"C",sep="_")
        colnames(drug_expr1) <- paste(colnames(drug_expr1),"T",sep="_")
        
        expr1 <- cbind(ctrl_1,drug_expr1)
        
        indiv <- factor(c(com_indiv,com_indiv))
        trt <- factor(c(rep("C",time=length(com_indiv)),c(rep("T",time=length(com_indiv)))),levels = c("C","T"))
        
        design1 <- model.matrix(~indiv + trt)
        fit2_1 <- lmFit(expr1,design1)
        fit2_1 <- eBayes(fit2_1)
        allDiff <- topTable(fit2_1,adjust = 'BH',coef = "trtT",n = Inf)
        all_fc_1 = as.data.frame(allDiff$logFC,colnames=drug,row.names = row.names(allDiff))
        write.table(all_fc_1,file=paste(output,drug,".txt",sep=""),sep="\t",quote=F)
        
    }
    }
  
  print(i)
  i=i+1
}


