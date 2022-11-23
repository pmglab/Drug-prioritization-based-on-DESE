# -*- coding: utf-8 -*-
# Input files:
# 1. ~/your_file_path/GSE92742_Broad_LINCS_inst_info.txt, instance information downloaded from GEO GSE92742;
# 2. ~/your_file_path/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx, level 3 data in LINCS downloaded from GEO GSE92742;
# 3. ~/your_file_path/GSE92742_Broad_LINCS_gene_info.txt, the information of the 12328 genes involved in LINCS downloaded from GEO GSE92742.

from sys import argv
script=argv

import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
from multiprocessing.dummy import Pool as ThreadPool

inst = pd.read_csv("~/your_file_path/GSE92742_Broad_LINCS_inst_info.txt",sep="\t")
ctrl = inst[inst.pert_type=="ctl_vehicle"]
ctrl_expr = parse("~/your_file_path/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",cid=ctrl.inst_id)
ctrl_expr_df = ctrl_expr.data_df
ctrl_expr_df["int_rid"] = ctrl_expr_df.index.map(int)

gene_info=pd.read_table("~/your_file_path/GSE92742_Broad_LINCS_gene_info.txt",sep="\t")
gene_ctrl = pd.merge(gene_info,ctrl_expr_df,left_on="pr_gene_id",right_on="int_rid")
del gene_ctrl["pr_gene_title"]
del gene_ctrl["pr_is_lm"]
del gene_ctrl["pr_is_bing"]
del gene_ctrl["int_rid"]
del gene_ctrl["pr_gene_id"]
gene_ctrl.set_index("pr_gene_symbol",inplace=True)
print("Converted rid to gene name!")
print(gene_ctrl.head())

id_list = ctrl.groupby(ctrl.cell_id)["inst_id"].apply(list)
aa =id_list.to_dict()
#print(aa)

def get_mean(drug):
    drug_mean = gene_ctrl[aa[drug]].mean(axis=1)
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

ctrl_expr = pd.DataFrame(drug_mean_dict,index=gene_ctrl.index)
ctrl_expr.to_csv("Mean_allLINCS.txt",sep="\t")
