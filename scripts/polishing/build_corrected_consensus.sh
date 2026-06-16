ml bedtools 
new_ver=$1

bedtools subtract -a v0.6_chrM.bed -b patch.bed > v0.6_include_for_polishing.bed
bcftools view -R v0.6_include_for_polishing.bed --no-version --threads 24 -Oz v0.6_snv_candidates.merfin-loose.vcf.gz > v0.6_snv_correction.vcf.gz
bcftools index v0.6_snv_correction.vcf.gz
bcftools concat --no-version --threads 24 -a -Oz v0.6_snv_correction.vcf.gz KOLF2.1Jv0.6.patch.vcf.gz  > v0.6_snv_sv_correction.vcf.gz
bcftools index v0.6_snv_sv_correction.vcf.gz
bcftools consensus -c v0.6_to_$new_ver.chain -HA -f v0.6_chrM.fa v0.6_snv_sv_correction.vcf.gz > $new_ver.fa
