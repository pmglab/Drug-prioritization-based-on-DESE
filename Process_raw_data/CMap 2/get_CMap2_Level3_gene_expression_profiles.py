from sys import argv
#'cell' denotes the cell line's name in CMap 2, such as NEU, NPC
#'time' denotes the perturbation time in CMap 2, such as 6, 24.
#'output' denotes the path to place the drug-induced gene expression profiles
script,cell,time,output=argv

import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress.parse import parse

#get approved or investigational drugs' names. Information can be accessed from DrugBank using r package 'dbparser'(https://cran.r-project.org/web/packages/dbparser/vignettes/dbparser.html)
drugbank_table = pd.read_table("./drugbank_drugName_group.txt",sep="\t")
appr_invst_drug = drugbank_table.name[(drugbank_table.group.str.contains("approved"))|(drugbank_table.group.str.contains("investigational"))].str.lower()

#load the instance information of GSE92742 (downloaded from GSE92742)
inst = pd.read_csv("./GSE92742_Broad_LINCS_inst_info.txt",sep="\t")

#extract informatin related to your 'cell'
inst_cell = inst[inst.cell_id==cell]

#extract information related to your 'cell' and 'time'
inst_cp = inst_cell[inst_cell.pert_iname.isin(appr_invst_drug)& (inst_cell.pert_time==int(time))]

#load instances' quality information from /GSE92742_Broad_LINCS_sig_metrics.txt(the distil_cc_q75 value).
qc_info = pd.read_table("./GSE92742_Broad_LINCS_sig_metrics.txt",sep="\t")

#extract the instance information with distil_cc_q75 > 0.15
qc_trt = qc_info[(qc_info.sig_id.str.contains("_"+cell+"_") & qc_info.sig_id.str.contains(str(time)+"H")) & (qc_info.pert_type == "trt_cp")]
qc_trt["plate"]=qc_info.sig_id.str.split("_",expand=True)[0]
qc_trt = qc_trt.sort_values(by=["plate","pert_iname","distil_cc_q75"],ascending=False)
qc_trt = qc_trt.drop_duplicates(subset=["plate","pert_iname"],keep="first")
print(qc_trt.iloc[:50,:])
qc_trt_high =qc_trt[qc_trt.distil_cc_q75>0.15]
print(qc_trt_high.iloc[:50,:])
qc_trt_high["ID"] = qc_trt_high.sig_id.str.split(":",expand=True)[0]+qc_trt_high.pert_id

inst_cp["ID"]=inst_cp.inst_id.str.split("_",expand=True)[0]+"_"+inst_cp.inst_id.str.split("_",expand=True)[1]+"_"+inst_cp.inst_id.str.split("_",expand=True)[2]+inst_cp.pert_id
high_qual_inst_cp = pd.merge(qc_trt_high,inst_cp,on="ID")

#extract the drug-induced gene expression profiles of your 'cell' and 'time'
cell_expr = parse("./GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",cid=high_qual_inst_cp.inst_id)
cell_expr_df = cell_expr.data_df
cell_expr_df["int_rid"] = cell_expr_df.index.map(int)

#extract the vehicle-induced gene expression profiles of your 'cell' and 'time' and also extract the instance information with high quality (distil_cc_q75 > 0.15)
ctrl = inst_cell[(inst_cell.pert_type=="ctl_vehicle") & (inst_cell.pert_time==int(time))]
ctrl["ID"]=ctrl.inst_id.str.split("_",expand=True)[0]+"_"+ctrl.inst_id.str.split("_",expand=True)[1]+"_"+ctrl.inst_id.str.split("_",expand=True)[2]+":"+ctrl.rna_well
qc_ctrl = qc_info[(qc_info.sig_id.str.contains("_"+cell+"_") & qc_info.sig_id.str.contains(str(time)+"H")) & (qc_info.pert_type == "ctl_vehicle")]
qc_ctrl["plate"]=qc_ctrl.sig_id.str.split("_",expand=True)[0]
qc_ctrl = qc_ctrl.sort_values(by=["plate","distil_cc_q75"],ascending=False)
qc_ctrl = qc_ctrl.drop_duplicates(subset=["plate","pert_iname"],keep="first")
qc_ctrl_high = qc_ctrl[qc_ctrl.distil_cc_q75>0.15]
ctrl_high = pd.merge(ctrl,qc_ctrl_high,left_on="ID",right_on="sig_id")

ctrl_expr = parse("./GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx",cid=ctrl_high.inst_id)
ctrl_expr_df = ctrl_expr.data_df
ctrl_expr_df["int_rid"] = ctrl_expr_df.index.map(int)

#convert gene rid in drug-iduced gene expression dataframe to gene name('GSE92742_Broad_LINCS_gene_info.txt' can be downloaded from GEO)
gene_info=pd.read_table("./GSE92742_Broad_LINCS_gene_info.txt",sep="\t")
gene = pd.merge(gene_info,cell_expr_df,left_on="pr_gene_id",right_on="int_rid")
del gene["pr_gene_title"]
del gene["pr_is_lm"]
del gene["pr_is_bing"]
del gene["int_rid"]
del gene["pr_gene_id"]
gene.set_index("pr_gene_symbol",inplace=True)
print("Converted rid to gene name!")

#convert gene rid in vehicle-iduced gene expression dataframe to gene name
gene_ctrl = pd.merge(gene_info,ctrl_expr_df,left_on="pr_gene_id",right_on="int_rid")
del gene_ctrl["pr_gene_title"]
del gene_ctrl["pr_is_lm"]
del gene_ctrl["pr_is_bing"]
del gene_ctrl["int_rid"]
del gene_ctrl["pr_gene_id"]
gene_ctrl.set_index("pr_gene_symbol",inplace=True)
print("Converted rid to gene name!")

#combine drug-induced and vehicle-induced gene expression profile dataframe
trt_ctrl = pd.merge(gene,gene_ctrl,left_index=True,right_index=True)
trt_ctrl.to_csv(output+cell+"_"+str(time)+"_pyschiatric_expr_high_quality_instance_distil_cc_q75_0.15.txt.txt",sep="\t")

#keep the instance information related to your 'cell' and 'time'
item = ["inst_id","pert_iname_x","pert_dose","rna_plate","pert_type_x"]
inst1 = pd.concat([high_qual_inst_cp[item],ctrl_high[item]],axis=0)
inst1.rename(columns={"pert_iname_x":"pert_iname","pert_type_x":"pert_type"},inplace=True)
inst1.to_csv(output+cell+"_"+str(time)+"psychiatric_inst_info_high_quality_instance_distil_cc_q75_0.15.txt",sep="\t")
