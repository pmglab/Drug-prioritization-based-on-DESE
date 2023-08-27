from sys import argv
script,cell,time=argv

import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress.parse import parse
from multiprocessing.dummy import Pool as ThreadPool
import time as T
start = T.clock()

#get approved or investigational drugs' names
drugbank = pd.read_table("resource/drugbank_drugName_group.txt",sep="\t")
app_drug = drugbank.name[(drugbank.group.str.contains("approved"))|(drugbank.group.str.contains("investigational"))].str.lower()

#download GSE92742_Broad_LINCS_inst_info.txt from GSE92742
inst = pd.read_csv("./GSE92742_Broad_LINCS_inst_info.txt",sep="\t")
inst_cell = inst[inst.cell_id==cell]
inst_cp = inst_cell[inst_cell.pert_iname.isin(app_drug)& (inst_cell.pert_time==int(time))]

#download GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx from GSE92742
cell_expr = parse("./GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",cid=inst_cp.inst_id)
cell_expr_df = cell_expr.data_df
cell_expr_df["int_rid"] = cell_expr_df.index.map(int)

ctrl = inst_cell[(inst_cell.pert_type=="ctl_vehicle") & (inst_cell.pert_time==int(time))]
ctrl_expr = parse("./GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",cid=ctrl.inst_id)
ctrl_expr_df = ctrl_expr.data_df
ctrl_expr_df["int_rid"] = ctrl_expr_df.index.map(int)

#convert drug expression file rid to gene name
#download GSE92742_Broad_LINCS_gene_info.txt from GSE92742
gene_info=pd.read_table("./GSE92742_Broad_LINCS_gene_info.txt",sep="\t")
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

#the full gene expression profiles of the user provided "cell" and "time".
trt_ctrl = pd.merge(gene,gene_ctrl,left_index=True,right_index=True)
trt_ctrl.to_csv("your_folder_path/"+cell+"_"+str(time)+"_expr.txt",sep="\t")

#the instance information of the user provided "cell" and "time".
inst1 = pd.concat([inst_cp,ctrl],axis=0)
inst1.to_csv("your_folder_path/"+cell+"_"+str(time)+"_inst_info.txt",sep="\t")