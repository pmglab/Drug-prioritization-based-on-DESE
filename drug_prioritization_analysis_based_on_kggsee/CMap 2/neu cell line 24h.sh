
#!/bin/bash
java -Xmx20g -jar kggsee_path/kggsee.jar --nt 10 --buildver hg19 --vcf-ref genotype_path/eur/1kg.phase3.v5.shapeit2.eur.hg19.chr_CHROM_.vcf.gz  --sum-file gwas_summary_path/SCZ/daner_PGC_SCZ_w3_90_0518d_eur.gz --filter-maf-le 0.05  --db-gene refgene  --gene-condi  --regions-out chr6:28477797-33448354  --expression-file ./drug_perturbation_profile_path/NEU_24H_drug_induced_foldchange_se.txt --out ./scz_lincs_24h_neu  --dese-permu-num 100  --multiple-testing bonf  --p-value-cutoff 0.05 --calc-selectivity false &&
java -Xmx20g -jar kggsee_path/kggsee.jar --nt 10 --buildver hg38 --vcf-ref genotype_path/eur_hg38/1kg.chr_CHROM_.shapeit2_integrated_snvindels_v2a_27022019.eur.GRCh38.phased.vcf.gz   --sum-file gwas_summary_path/Depression/PGC_UKB_depression_genome-wide_update_hg38.txt.result.txt.gz --filter-maf-le 0.05  --db-gene refgene  --gene-condi  --regions-out chr6:28510120-33480577  --expression-file ./drug_perturbation_profile_path/NEU_24H_drug_induced_foldchange_se.txt --out ./mdd_lincs_24h_neu --calc-selectivity false --dese-permu-num 100  --multiple-testing benfdr  --p-value-cutoff 0.05 &&
java -Xmx20g -jar kggsee_path/kggsee.jar --chrom-col CHROM --pos-col POS --p-col PVAL --nt 10 --vcf-ref genotype_path/eur/1kg.phase3.v5.shapeit2.eur.hg19.chr_CHROM_.vcf.gz   --sum-file gwas_summary_path/BD/pgc-bip2021-BDI_delete_#.vcf.tsv --filter-maf-le 0.05  --db-gene refgene  --gene-condi   --regions-out chr6:28477797-33448354 --expression-file ./drug_perturbation_profile_path/NEU_24H_drug_induced_foldchange_se.txt --out ./bdi_lincs_24h_neu --calc-selectivity false  --dese-permu-num 100  --multiple-testing benfdr  --p-value-cutoff 0.01 
