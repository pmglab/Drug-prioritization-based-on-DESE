library(GEOquery)

#GSE119291_family.soft.gz was downloaded from GSE119291
data <- getGEO(filename="./GSE119291_family.soft.gz",GSEMatrix=F)
sample <- data@header$sample_id

#extact data from GEO
expr_all = data@gsms[[sample[1]]]@dataTable@table
colnames(expr_all)[2] <- sample[1]

for (i in 2:5674){
    expr <- data@gsms[[sample[i]]]@dataTable@table
    if (dim(expr)[1]>0){
        if (sum(data@gsms[[sample[1]]]@dataTable@table$ID_REF == expr$ID_REF) == 22268){
        colnames(expr)[2] <- sample[i]
        expr_all <- cbind(expr_all,expr)
        expr_all <- expr_all[,-(i+1)]
        }
    }else{print(sample[i])}
}
write.table(expr_all,file="./GSE19291_all_expr.txt",sep="\t",row.names=F,quote=F)

expr_all = read.table("./GSE19291_all_expr.txt",sep="\t",header=T)

#add gene symbol info (GPL25480	Genometry L1000™ Expression Profiling)
gene <- read.table("./gene_annotation.txt",sep="\t",header=F)
expr_gene <- merge(expr_all,gene,by.x="ID_REF",by.y="V1")

#merge same gene expression
expr_mean <- aggregate(expr_gene[,2:5675],by=list(expr_gene$V2),FUN=mean)
write.table(expr_mean,file="./GSE19291_all_expr_geneSymbol.txt",sep="\t",quote=F)

#get all normal hipsc npc and sz npc gene expression
#get from GEO website(GSE119291)
sample_info = read.table("./sample_info.txt",sep="\t",header=T)
print(head(sample_info))

norm_info = sample_info[which(sample_info$Source.name == "hiPSC_control"),]$Accession
sz_info = sample_info[which(sample_info$Source.name == "hiPSC_SZ"),]$Accession

norm_expr = expr_mean[,norm_info]
sz_expr = expr_mean[,sz_info]

write.table(norm_expr,file="./all_hips_ctl_expression.txt",sep="\t",row.names=T,quote=F)
write.table(sz_expr,file="./all_hips_sz_expression.txt",sep="\t",row.names=T,quote=F)
