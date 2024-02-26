
args <- commandArgs(trailingOnly = TRUE)
raw_expr_path <- args[1] #Get this information from GEO website(GSE119291),i.e.,GSE119291_family.soft.gz
sample_info_path <- args[2] #Get https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE119291&platform=GPL25480
hipsc_type <- args[3]#"hiPSC_control" or "hiPSC_SZ"
output <- args[4]#Folders to place result files

library(limma)
library(stringr)

#load the gene expression profiles
data <- read.table(raw_expr_path,header=T,sep="\t",row.names = 1,check.names=FALSE)

#load the sample information of hiPSC_control or hiPSC_SZ
pheno = read.table(sample_info_path,sep="\t",header=T)
pheno = pheno[which(pheno$Source.name == hipsc_type),]

#remove the two non-European individuals' cell lines
pheno = pheno[(which(pheno$"Cell.id"!="3084-2-2")),]
pheno = pheno[(which(pheno$"Cell.id"!="3234-2-4")),]
pheno["drug_dose"] = paste(pheno$Perturbagen,pheno$Dosage,sep="_")
pheno["new_batch"] = paste(pheno$Cell.id,pheno$Batch,sep="_")

#statistic the number of perturbagens in total
all_drugs = unique(pheno$drug_dose)
print(paste("There are",length(all_drugs),"drug-dose pair",sep=" "))
print(all_drugs)

#control vehicle-induced gene expression profiles
ctrl_expr = data[,pheno[which(pheno$Perturbation.type=="vehicle"),]$Accession]#raw control vehicle-induced gene expression profiles
ctrl_expr_pheno = merge(pheno,as.data.frame(t(ctrl_expr)),by.x="Accession",by.y=0)# add sample information for control expression
ctrl_expr_pheno = ctrl_expr_pheno[,c("new_batch",colnames(t(ctrl_expr)))]# keep individual information for control expression
ctrl_expr_pheno = aggregate(ctrl_expr_pheno[,2:dim(ctrl_expr_pheno)[2]],by=list(ctrl_expr_pheno$new_batch),FUN = mean)# merge control expression from same individual
ctrl_expr_pheno = as.data.frame(ctrl_expr_pheno[,2:dim(ctrl_expr_pheno)[2]],row.names = as.character(ctrl_expr_pheno$Group.1))#individual id as the index column
ctrl_merged = t(ctrl_expr_pheno)

#A function to compute the standard deviation of the log2foldchanges of a certain gene induced by a treatment
std <- function(x){return (sd(x)/sqrt(length(x)))}

#gene differential expression analysis based on limma
i=1
for (dg in all_drugs[1:length(all_drugs)]){
  if (!str_detect(dg,"DMSO")){#iterate all drugs
    print(dg)
    drug_expr = data[,as.character(pheno$Accession[pheno$drug_dose==dg])]
    
    if (class(drug_expr)!= "numeric"){
        com_indiv = intersect(colnames(ctrl_merged),unique(as.character(pheno$new_batch[pheno$drug_dose==dg])))#select individuals treated with vehicle and drugs
        
        #get drug-induced gene expression profiles
        drug = merge(pheno,t(drug_expr),by.x="Accession",by.y=0)
        drug = drug[,c("new_batch",colnames(t(drug_expr)))]
        drug = aggregate(drug[,2:dim(drug)[2]],by=list(drug$new_batch),FUN = mean)
        drug = as.data.frame(drug[,2:dim(drug)[2]],row.names = as.character(drug$Group.1))
        drug_expr = t(drug)
        
        ctrl_1 = ctrl_merged[,com_indiv]
        drug_expr1 = drug_expr[,com_indiv]
        
        fc = drug_expr1 - ctrl_1
        se = apply(fc,MARGIN = 1,FUN = std)
        se = as.data.frame(se)
        #print(head(se))

        #get drug-induced and vehicle-induced gene expression profiles from same individuals
        ctrl_expr = ctrl_merged[,com_indiv]
        drug_expr = drug_expr[,com_indiv]
        
        #design control and treatment group 
        colnames(ctrl_expr) <- paste(colnames(ctrl_expr),"C",sep="_")
        colnames(drug_expr) <- paste(colnames(drug_expr),"T",sep="_")
        
        expr <- merge(ctrl_expr,drug_expr,by="row.names")
        expr <- as.data.frame(expr[,2:dim(expr)[2]],row.names=expr$Row.names)
        print(head(expr)) 
        indiv <- factor(c(com_indiv,com_indiv))
        trt <- factor(c(rep("C",time=length(com_indiv)),c(rep("T",time=length(com_indiv)))),levels = c("C","T"))

        #perform differential expression analysis
        design <- model.matrix(~indiv + trt)
        fit <- lmFit(expr,design)
        fit <- eBayes(fit)
        allDiff <- topTable(fit,adjust = 'BH',coef = "trtT",n = Inf)

        #combine the log2foldchanges and se of each gene
        allDiff <- merge(se,allDiff,by.x=0,by.y=0)
        all_fc_1 = as.data.frame(allDiff[,c("Row.names","logFC","se")])
        colnames(all_fc_1)[2]=dg
        colnames(all_fc_1)[3]=paste(dg,".SE",sep="")

        #write the result of log2foldchange and se to a file
        write.table(all_fc_1,file=paste(output,dg,".txt",sep=""),sep="\t",quote=F,row.names=F)
        
    }
    }
  
  print(i)
  i=i+1
}


