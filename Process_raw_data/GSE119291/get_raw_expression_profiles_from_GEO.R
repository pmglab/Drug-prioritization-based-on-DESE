
args <- commandArgs(trailingOnly = TRUE)
family_data_path <- args[1] #Get this information from GEO website(GSE119291),i.e.,GSE119291_family.soft.gz
sample_info_path <- args[2] #Get https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE119291&platform=GPL25480
gene_annotation_file_path <- args[3]#Affymetrix (HgU133A) probe information
output <- args[4] #Folders to place result files

library(GEOquery)
data <- getGEO(filename=family_data_path,GSEMatrix=F)
sample <- data@header$sample_id

#extact all samples' data from GEO family data
expr_all = data@gsms[[sample[1]]]@dataTable@table
colnames(expr_all)[2] <- sample[1]
#5674 drug-induced gene expression profiles generated on GPL25480 platform were used in this study
for (i in 2:5674){
    expr <- data@gsms[[sample[i]]]@dataTable@table
    if (dim(expr)[1]>0){
        if (sum(data@gsms[[sample[1]]]@dataTable@table$ID_REF == expr$ID_REF) == 22268){
        colnames(expr)[2] <- sample[i]
        expr_all <- cbind(expr_all,expr)
        expr_all <- expr_all[,-(i+1)]# delete the redundant ID columns
        }
    }else{print(sample[i])}
}

#covert affymetrix id to correspondig gene symbol
#Affymetrix (HgU133A) probe information was downloaded
gene <- read.table(gene_annotation_file_path,sep="\t",header=F)
expr_gene <- merge(expr_all,gene,by.x="ID_REF",by.y="V1")

#merge the expression values of multiple probes corresponding to the same gene
expr_mean <- aggregate(expr_gene[,2:5675],by=list(expr_gene$V2),FUN=mean)
rownames(expr_mean)<-expr_mean[,1]
expr_mean <- expr_mean[,-1]

#get sample information to further extract normal hipsc and sz gene expression profiles, respectively
sample_info = read.table(sample_info_path,sep="\t",header=T)#get https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE119291&platform=GPL25480

#extract the Accession ID of Control NPC line and hiPSC NPC line's drug- and control-vehicle-induced gene expression profiles
ctrl_npc_info = sample_info[which(sample_info$Source.name == "hiPSC_control"),]$Accession
sz_npc_info = sample_info[which(sample_info$Source.name == "hiPSC_SZ"),]$Accession

#extract the Control NPC line and hiPSC NPC line's drug- and control-vehicle-induced gene expression profiles
ctrl_npc_expr = expr_mean[,ctrl_npc_info]
sz_npc_expr = expr_mean[,sz_npc_info]

#write the drug- and control-vehicle- induced gene expression profiles to files.
write.table(ctrl_npc_expr,file=paste(output,"gse119291_ctrl_hipsc_expression.txt",sep=""),sep="\t",row.names=T,quote=F)
write.table(sz_npc_expr,file=paste(output,"gse119291_sz_hipsc_expression.txt",sep=""),sep="\t",row.names=T,quote=F)
