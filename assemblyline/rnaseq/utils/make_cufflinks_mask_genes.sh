awk '$1 == "MT"' genes.gtf > cufflinks_mask_genes.gtf
awk '$2 == "Mt_rRNA"' genes.gtf >> cufflinks_mask_genes.gtf
awk '$2 == "Mt_tRNA"' genes.gtf >> cufflinks_mask_genes.gtf
awk '$2 == "Mt_tRNA_pseudogene"' genes.gtf >> cufflinks_mask_genes.gtf
awk '$2 == "rRNA"' genes.gtf >> cufflinks_mask_genes.gtf
awk '$2 == "rRNA_pseudogene"' genes.gtf >> cufflinks_mask_genes.gtf
cat cufflinks_mask_genes.gtf | sort | uniq > a
mv a cufflinks_mask_genes.gtf
