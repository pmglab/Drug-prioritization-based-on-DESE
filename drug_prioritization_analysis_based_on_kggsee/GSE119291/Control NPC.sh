
#!/bin/bash
java -Xmx20g -jar kggsee_path/kggsee.jar --nt 10 --buildver hg19 --vcf-ref genotype_path/gty/eur/1kg.phase3.v5.shapeit2.eur.hg19.chr_CHROM_.vcf.gz --sum-file summary_file_path/gwas_file/SCZ/daner_PGC_SCZ_w3_90_0518d_eur.gz --filter-maf-le 0.05  --db-gene refgene  --gene-condi  --regions-out chr6:28477797-33448354  --expression-file ./expr_prof_for_kggsee/hipsc_ctrl_with_se_drug_induced_foldchange.txt --out /scz_geo_ctrl --dese-permu-num 100  --multiple-testing bonf  --p-value-cutoff 0.05 --calc-selectivity false &&
java -Xmx20g -jar kggsee_path/kggsee.jar --nt 10 --buildver hg38 --vcf-ref genotype_path/gty/eur_hg38/1kg.chr_CHROM_.shapeit2_integrated_snvindels_v2a_27022019.eur.GRCh38.phased.vcf.gz  --sum-file summary_file_path/gwas_file/Depression/PGC_UKB_depression_genome-wide_update_hg38.txt.result.txt.gz --filter-maf-le 0.05  --db-gene refgene  --gene-condi  --regions-out chr6:28510120-33480577  --expression-file ./expr_prof_for_kggsee/hipsc_ctrl_with_se_drug_induced_foldchange.txt --out /mdd_geo_ctrl  --dese-permu-num 100  --multiple-testing benfdr  --p-value-cutoff 0.05 --calc-selectivity false &&
java -Xmx20g -jar kggsee_path/kggsee.jar --chrom-col CHROM --pos-col POS --p-col PVAL --nt 10 --vcf-ref genotype_path/gty/eur/1kg.phase3.v5.shapeit2.eur.hg19.chr_CHROM_.vcf.gz --sum-file summary_file_path/gwas_file/BD/pgc-bip2021-BDI_delete_#.vcf.tsv --filter-maf-le 0.05  --db-gene refgene  --gene-condi   --regions-out chr6:28477797-33448354  --expression-file ./expr_prof_for_kggsee/hipsc_ctrl_with_se_drug_induced_foldchange.txt --out /bdi_geo_ctrl --dese-permu-num 100  --multiple-testing benfdr  --p-value-cutoff 0.01 --calc-selectivity false
