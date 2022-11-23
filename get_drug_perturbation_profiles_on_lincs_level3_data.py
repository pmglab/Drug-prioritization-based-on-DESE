# 1. celltype: on of the result files of DESE. 
#    This file contains three colums, in which the first is cell line name and the second is the association p-value of the cell line and the investigated disorder.
# 2. dise: a str, representing the name of the investigated disorder.

# This script also contains four other input files:
# 1. ~/your_file_path/drugbank_drugName_group.txt, a file containing two columns, the first is drug name, the second is the "Group" of the drug in DrugBank, such as "investigational" and "approved".
# 2. ~/your_file_path/GSE92742_Broad_LINCS_inst_info.txt, instance information downloaded from GEO GSE92742;
# 3. ~/your_file_path/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx, level 3 data in LINCS downloaded from GEO GSE92742;
# 4. ~/your_file_path/GSE92742_Broad_LINCS_gene_info.txt, the information of the 12328 genes involved in LINCS downloaded from GEO GSE92742.

from sys import argv
script,celltype,dise=argv

import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress.parse import parse
from multiprocessing.dummy import Pool as ThreadPool
import time as T
start = T.clock()

#approved drugs
drugbank = pd.read_table("~/your_file_path/drugbank_drugName_group.txt",sep="\t")
app_drug = drugbank.name[(drugbank.group.str.contains("approved"))|(drugbank.group.str.contains("investigational"))].str.lower()

cell = pd.read_table(celltype,sep="\t")
sigcell = cell.TissueName[cell.RobustRegressionZ<(0.05/cell.shape[0])]
print(sigcell)

inst = pd.read_csv("~/your_file_path/GSE92742_Broad_LINCS_inst_info.txt",sep="\t")
inst_ctl = inst[inst.cell_id.isin(sigcell)]
inst_cp = inst_ctl[inst_ctl.pert_iname.isin(app_drug)]

cell_expr = parse("~/your_file_path/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",cid=inst_cp.inst_id)
cell_expr_df = cell_expr.data_df
cell_expr_df["int_rid"] = cell_expr_df.index.map(int)
#print(cell_expr_df.head())

ctrl = inst_ctl[inst_ctl.pert_type=="ctl_vehicle"]
ctrl_expr = parse("~/your_file_path/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",cid=ctrl.inst_id)
ctrl_expr_df = ctrl_expr.data_df
ctrl_expr_df["int_rid"] = ctrl_expr_df.index.map(int)

#convert drug expression file rid to gene name
gene_info=pd.read_table("~/your_file_path/GSE92742_Broad_LINCS_gene_info.txt",sep="\t")
gene = pd.merge(gene_info,cell_expr_df,left_on="pr_gene_id",right_on="int_rid")
del gene["pr_gene_title"]
del gene["pr_is_lm"]
del gene["pr_is_bing"]
del gene["int_rid"]
del gene["pr_gene_id"]
gene.set_index("pr_gene_symbol",inplace=True)
print("Converted rid to gene name!")

#convert control expression file rid to gene name
gene_ctrl = pd.merge(gene_info,ctrl_expr_df,left_on="pr_gene_id",right_on="int_rid")
del gene_ctrl["pr_gene_title"]
del gene_ctrl["pr_is_lm"]
del gene_ctrl["pr_is_bing"]
del gene_ctrl["int_rid"]
del gene_ctrl["pr_gene_id"]
gene_ctrl.set_index("pr_gene_symbol",inplace=True)
print("Converted rid to gene name!")

#get log2 FC
gene_ctrl["Mean"]=gene_ctrl.mean(axis=1)
print(gene_ctrl.shape)

id_list = inst_cp.groupby(inst_cp.pert_iname)["inst_id"].apply(list)
aa =id_list.to_dict()
print(aa.keys())
print("There are "+str(len(aa.keys()))+" drugs.")

raw_mean = list()
def get_mean(drug):
    drug_mean = gene[aa[drug]].mean(axis=1)
    drug_mean = drug_mean.tolist()
    drug_mean.insert(0,drug)
    return(drug_mean)

pool = ThreadPool(10)
raw_mean = pool.map(get_mean, aa.keys())
pool.close()
pool.join()

drug_mean_dict = dict()
for item in raw_mean:
    drug_mean_dict[item[0]]=item[1:]

fc_expr2 = pd.DataFrame(drug_mean_dict,index=gene.index)

logGene = (fc_expr2+1).applymap(np.log2)
fc_expr1 = logGene.sub(np.log2(gene_ctrl.Mean+1),axis=0)
fc_expr2 = -fc_expr1

fc_expr1.to_csv("allTime_mean_"+dise+"_signif_cell_line's_trt_raw_level3_perturbation_profile.txt",sep="\t")
fc_expr2.to_csv("allTime_mean_"+dise+"_signif_cell_line's_trt_reverse_level3_perturbation_profile.txt",sep="\t")

elapsed = (T.clock() - start)
print("Time used:",elapsed)

