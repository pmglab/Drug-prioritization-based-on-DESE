This page provided the main code, basic data and scripts used in our drug prioritization analysis based on DESE(Genome Biol 20, 233 (2019). https://doi.org/10.1186/s13059-019-1801-5).

### **Step 1: Process the raw data and get the drug- and control-vehicle -induced gene expression profiles**
1)./Process_raw_data/GSE119291/get_raw_expression_profiles_from_GEO.R was used to get the raw drug-induced gene expression profiles of Control NPC and SZ NPC line from GSE119291.

2)./Process_raw_data/CMap2/get_CMap2_Level3_gene_expression_profiles.py was used to get the raw drug-induced gene expression profiles of certain cell line at certain time from CMap 2.

### **Step 2: Compute the genes' log2foldchanges induced by drugs (drug perturbation profiles)**
1)./get_drug_perturbation_profiles/GSE119291/Differential_gene_expression_analysis.R was used to perform the differential expression analysis based on the drug-induced gene expression profiles of Control NPC and SZ NPC line. The result was put into KGGSEE to perform selective drug perturbation analysis.

2)./get_drug_perturbation_profiles/CMap 2/Differential_gene_expression_analysis.R was used to perform the differential expression analysis based on the drug-induced gene expression profiles of certain cell line in CMap 2. The result was put into KGGSEE to perform selective drug perturbation analysis.

### **Step 3: Prioritize the candidate drugs for common psychiatric disorders based on DESE**
The shell scripts under the "drug_prioritization_analysis_based_on_kggsee" folder were used to perform the drug prioritization analysis based on DESE (implemented in our software platform KGGSEE, https://pmglab.top/kggsee/#/). A detailed description of the parameters of options can be seen at https://kggsee.readthedocs.io/en/latest/detailed_document.html#dese-for-drug-prioritization-analysis.

Note: The parameter of option "--vcf-ref" in shell scripts under the "drug_prioritization_analysis_based_on_kggsee" folder was the VCF files of genotypes sampled from EUR/EAS/AFR panel from the 1000 Genomes Project (phase 3) and can be downloaded from http://pmglab.top/genotypes/#/. The usage of "--vcf-ref" is "--vcf-ref 1kg.phase3.v5.shapeit2.eur.hg19.chr_CHROM_.vcf.gz".
