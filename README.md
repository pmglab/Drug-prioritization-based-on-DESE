The main code, basic data and script used in our drug repositioning research based on DESE.

./Process_raw_data/get_GSE119291_expr.R was used to get the raw gene expression profiles of Control NPC and SZ NPC from GSE119291.

./Process_raw_data/get_LINCS_level3_gene_expr.R was used to get the raw gene expression profiles of certain cell line at certain time from LINCS.

./differential_expression_analysis/GSE119291.R was used to perform the differential expression analysis based on the drug-induced gene expression profiles of Control NPC and SZ NPC. The result was put into KGGSEE to perform selective drug perturbation analysis.

./differential_expression_analysis/lincs.R was used to perform the differential expression analysis based on the drug-induced gene expression profiles of certain cell line in LINCS. The result was put into KGGSEE to perform selective drug perturbation analysis.

Note: the parameter of option "--vcf-ref" in "kggsee_for_drug_repositioning_research.sh" was the VCF files of genotypes sampled from EUR panel from the 1000 Genomes Project (phase 3) and can be downloaeded from http://pmglab.top/genotypes/#/. The usage of "--vcf-ref" is "--vcf-ref 1kg.phase3.v5.shapeit2.eur.hg19.chr_CHROM_.vcf.gz".
