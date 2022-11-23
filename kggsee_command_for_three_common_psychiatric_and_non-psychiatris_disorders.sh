java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar  --nt 12  --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/BD/pgc-bip2021-BDI_delete_#.vcf.tsv --filter-maf-le 0.05   --gene-finemapping --regions-out chr6:27477797-34448354 --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_BDI_fdr_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/bdi_raw_all_5kb_approved_drugs   --multiple-testing benfdr --p-value-cutoff 0.05 --chrom-col CHROM --p-col PVAL --pos-col POS --calc-selectivity --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar  --nt 12 --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/SCZ/daner_PGC_SCZ_w3_90_0518d_eur.gz --filter-maf-le 0.05   --gene-finemapping   --regions-out chr6:27477797-34448354 --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/scz_raw_all_5kb_approved_drugs_bonf0.05 --multiple-testing bonf --p-value-cutoff 0.05 --expression-file ~/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt --calc-selectivity --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar  --nt 12  --saved-ref /home/lxy/02ECS_eqtl/result/KGGSEE_maf05/VCFRefhg38 --buildver hg38 --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/Depression/PGC_UKB_depression_genome-wide_update_hg38.txt.result.txt.gz --buildver hg38 --filter-maf-le 0.05 --gene-finemapping  --regions-out chr6:27477797-34448354 --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/depression_raw_all_5kb_cell_approved_benfdr0.05 --multiple-testing benfdr --p-value-cutoff 0.05 --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_MDD2019_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt --calc-selectivity --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 12   --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/BD/pgc-bip2021-BDI_delete_#.vcf.tsv --filter-maf-le 0.05   --gene-finemapping --regions-out chr6:27477797-34448354 --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_BDI_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt  --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/bdi_reverse_all_5kb_approved_drugs   --multiple-testing benfdr --p-value-cutoff 0.05 --chrom-col CHROM --p-col PVAL --pos-col POS --calc-selectivity --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 12  --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/SCZ/daner_PGC_SCZ_w3_90_0518d_eur.gz --filter-maf-le 0.05   --gene-finemapping   --regions-out chr6:27477797-34448354 --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/scz_reverse_all_5kb_approved_drugs_bonf0.05 --multiple-testing bonf  --p-value-cutoff 0.05 --expression-file ~/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt --calc-selectivity --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 12   --saved-ref /home/lxy/02ECS_eqtl/result/KGGSEE_maf05/VCFRefhg38 --buildver hg38 --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/Depression/PGC_UKB_depression_genome-wide_update_hg38.txt.result.txt.gz --buildver hg38 --filter-maf-le 0.05 --gene-finemapping  --regions-out chr6:27477797-34448354 --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/depression_reverse_all_5kb_cell_approved_benfdr0.05 --multiple-testing benfdr --p-value-cutoff 0.05 --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_MDD2019_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt --calc-selectivity --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar  --nt 12 --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/SCZ/daner_PGC_SCZ_w3_90_0518d_eur.gz --filter-maf-le 0.05   --gene-finemapping   --regions-out chr6:27477797-34448354 --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/scz_raw_all_5kb_approved_drugs_bonf0.05_BF_eQTL --multiple-testing bonf --p-value-cutoff 0.05 --expression-file ~/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt --calc-selectivity --eqtl-file /home/lxy/02ECS_eqtl/eQTLs/01_gtexv8_eQTL/02_eur_eQTL/EUR_gene_eqtl_hg19/Brain_Frontal_Cortex_BA9_eur_v8_tmm_p01.gene.hg19.cov.eqtl.txt.gz --filter-eqtl-p 0.01 --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 12  --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/SCZ/daner_PGC_SCZ_w3_90_0518d_eur.gz --filter-maf-le 0.05   --gene-finemapping   --regions-out chr6:27477797-34448354 --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/scz_reverse_all_5kb_approved_drugs_bonf0.05_BF_eQTL --multiple-testing bonf  --p-value-cutoff 0.05 --expression-file ~/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt --calc-selectivity --eqtl-file /home/lxy/02ECS_eqtl/eQTLs/01_gtexv8_eQTL/02_eur_eQTL/EUR_gene_eqtl_hg19/Brain_Frontal_Cortex_BA9_eur_v8_tmm_p01.gene.hg19.cov.eqtl.txt.gz --filter-eqtl-p 0.01 --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar  --nt 12  --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/BD/pgc-bip2021-BDI_delete_#.vcf.tsv --filter-maf-le 0.05   --gene-finemapping --regions-out chr6:27477797-34448354 --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_BDI_fdr_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/bdi_raw_all_5kb_approved_drugs_BF_eQTL --multiple-testing benfdr --p-value-cutoff 0.05 --chrom-col CHROM --p-col PVAL --pos-col POS --calc-selectivity --eqtl-file /home/lxy/02ECS_eqtl/eQTLs/01_gtexv8_eQTL/02_eur_eQTL/EUR_gene_eqtl_hg19/Brain_Frontal_Cortex_BA9_eur_v8_tmm_p01.gene.hg19.cov.eqtl.txt.gz --filter-eqtl-p 0.01 --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 12   --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/BD/pgc-bip2021-BDI_delete_#.vcf.tsv --filter-maf-le 0.05   --gene-finemapping --regions-out chr6:27477797-34448354 --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_BDI_fdr_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt  --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/bdi_reverse_all_5kb_approved_drugs_BF_eQTL  --multiple-testing benfdr --p-value-cutoff 0.05 --chrom-col CHROM --p-col PVAL --pos-col POS --calc-selectivity --eqtl-file /home/lxy/02ECS_eqtl/eQTLs/01_gtexv8_eQTL/02_eur_eQTL/EUR_gene_eqtl_hg19/Brain_Frontal_Cortex_BA9_eur_v8_tmm_p01.gene.hg19.cov.eqtl.txt.gz --filter-eqtl-p 0.01 --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar  --nt 12  --saved-ref /home/lxy/02ECS_eqtl/result/KGGSEE_maf05/VCFRefhg38 --buildver hg38 --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/Depression/PGC_UKB_depression_genome-wide_update_hg38.txt.result.txt.gz --buildver hg38 --filter-maf-le 0.05 --gene-finemapping  --regions-out chr6:27477797-34448354 --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/depression_raw_all_5kb_cell_approved_benfdr0.05_BAN_eQTL --multiple-testing benfdr --p-value-cutoff 0.05 --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_MDD2019_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt --calc-selectivity --eqtl-file /home/lxy/02ECS_eqtl/eQTLs/01_gtexv8_eQTL/02_eur_eQTL/EUR_gene_eqtl_hg38/Brain_Anterior_cingulate_cortex_BA24_eur_v8_tmm_p01.gene.hg38.cov.eqtl.txt.gz --filter-eqtl-p 0.01 --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 12   --saved-ref /home/lxy/02ECS_eqtl/result/KGGSEE_maf05/VCFRefhg38 --buildver hg38 --sum-file /home/lxy/02ECS_eqtl/result/gwas_file/Depression/PGC_UKB_depression_genome-wide_update_hg38.txt.result.txt.gz --buildver hg38 --filter-maf-le 0.05 --gene-finemapping  --regions-out chr6:27477797-34448354 --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/depression_reverse_all_5kb_cell_approved_benfdr0.05_BAN_eQTL --multiple-testing benfdr --p-value-cutoff 0.05 --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_MDD2019_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt --calc-selectivity --eqtl-file /home/lxy/02ECS_eqtl/eQTLs/01_gtexv8_eQTL/02_eur_eQTL/EUR_gene_eqtl_hg38/Brain_Anterior_cingulate_cortex_BA24_eur_v8_tmm_p01.gene.hg38.cov.eqtl.txt.gz --filter-eqtl-p 0.01 --only-hgnc-gene &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 10 --pfile /home/lxy/02ECS_eqtl/result/gwas_file/T2D/Xue_et_al_T2D_META_Nat_Commun_2018.gz --chrom-col CHR --pos-col BP --p-col P --gene-finemapping --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/T2D_scz_signifCell_allTime_level3noabs_approvedInvesti_add_one --filter-maf-le 0.05 --only-hgnc-gene --p-value-cutoff 0.05 --multiple-testing bonf --regions-out chr6:27477797-34448354 --calc-selectivity --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 10 --pfile /home/lxy/02ECS_eqtl/result/gwas_file/T2D/Xue_et_al_T2D_META_Nat_Commun_2018.gz --chrom-col CHR --pos-col BP --p-col P --gene-finemapping --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/T2D_scz_signifCell_allTime_level3reverse_approvedInvesti_add_one --filter-maf-le 0.05 --only-hgnc-gene --p-value-cutoff 0.05 --multiple-testing bonf --regions-out chr6:27477797-34448354 --calc-selectivity --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 10 --pfile /home/lxy/02ECS_eqtl/result/gwas_file/RA/RA_GWASmeta_European_v2.txt --chrom-col Chr --pos-col Position --p-col P-val --gene-finemapping --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --out DR_psy/replication_test/RA_scz_signifCell_allTime_level3noabs_approvedInvesti_add_one --filter-maf-le 0.05 --only-hgnc-gene --p-value-cutoff 0.05 --multiple-testing bonf --regions-out chr6:27477797-34448354 --calc-selectivity --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 10 --pfile /home/lxy/02ECS_eqtl/result/gwas_file/RA/RA_GWASmeta_European_v2.txt --chrom-col Chr --pos-col Position --p-col P-val --gene-finemapping --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/RA_scz_signifCell_allTime_level3reverse_approvedInvesti_add_one --filter-maf-le 0.05 --only-hgnc-gene --p-value-cutoff 0.05 --multiple-testing bonf --regions-out chr6:27477797-34448354 --calc-selectivity --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt && 
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 10 --pfile /home/lxy/02ECS_eqtl/result/gwas_file/CAD/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz --chrom-col chr --pos-col bp_hg19 --p-col p-value_gc --gene-finemapping --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/CAD_scz_signifCell_allTime_level3no_abs_approvedInvesti_add_one --filter-maf-le 0.05 --only-hgnc-gene --p-value-cutoff 0.05 --multiple-testing bonf --regions-out chr6:27477797-34448354 --calc-selectivity --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_no_abs_level3_appInvesti_add_one.txt &&
java8 -Xmx10g -jar ~/kggsee/kggsee_old.jar --nt 10 --pfile /home/lxy/02ECS_eqtl/result/gwas_file/CAD/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz --chrom-col chr --pos-col bp_hg19 --p-col p-value_gc --gene-finemapping --saved-ref /home/lxy/02ECS_eqtl/result/isoformQTLs/tissue_selective/SCZ/SCZ_gene --out /home/lxy/03Drug_Reposition/DR_psy/replication_test/CAD_scz_signifCell_allTime_level3reverse_approvedInvesti_add_one --filter-maf-le 0.05 --only-hgnc-gene --p-value-cutoff 0.05 --multiple-testing bonf --regions-out chr6:27477797-34448354 --calc-selectivity --expression-file /home/lxy/03Drug_Reposition/DR_psy/expr_foldChange/parallel_Mean_allTime_SCZ_signif_cell_line_trt_reverse_level3_appInvesti_add_one.txt 
