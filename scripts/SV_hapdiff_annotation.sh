# Processing hapdiff variants
### Downloaded gap.bed and segdup.bed from GIAB
### Kept only the centromeres from the gap.bed
### Made the centromere and segdup beds into bed3
### Combined both files to have an exclude region

## 1. Filter hapdiff vcf to keep only SVs minimum 50bp
bcftools view \
  -i 'abs(INFO/SVLEN)>=50' \
  hapdiff_phased.vcf.gz \
  -Oz -o hapdiff_phased.50bp.vcf.gz

## 2. Filter out only centromeres
bcftools view \
  --targets-file ^centromere.bed3 \
  hapdiff_phased.50bp.vcf.gz \
  -Oz \
  -o hapdiff_phased.50bp.noCentromere.vcf.gz

bcftools index hapdiff_phased.50bp.noExcluded.vcf.gz

## 3. Annotate with AnnotSV
conda activate annotsv
$PATH/AnnotSV \
  -SVinputFile hapdiff_phased.50bp.noCentromeres.vcf.gz \
  -genomeBuild GRCh38 \
  -outputFile hapdiff.50bp.nocentromere.tsv


## 4. Try to summarize, removing SV sequence
awk -F'\t' -v COL=10 '
NR==1{
    for(i=1;i<=NF;i++)
        if(i!=COL)
            printf "%s%s",$i,(i==NF?"\n":"\t")
    next
}
{
key=$1 FS $2 FS $3 FS $4
genes[key]=(genes[key]?genes[key]","$5:$5)
if(!seen[key]++){
    row[key]=$0
}
}
END{
for(k in row){
    split(row[k],a,"\t")
    a[5]=genes[k]
    for(i=1;i<=length(a);i++)
        if(i!=COL)
            printf "%s%s",a[i],(i==length(a)?"\n":"\t")
}
}' hapdiff.50bp.nocentromere.annotated.tsv > hapdiff.50bp.nocentromere.summary.tsv


## Annovar
tar -xzf annovar.latest.tar.gz

convert2annovar.pl -format vcf4 hapdiff_phased.50bp.noCentromere.vcf.gz > hapdiff_phased.50bp.noCentromere.avinput
cd annovar
table_annovar.pl ../hapdiff_phased.50bp.noCentromere.avinput humandb/ \
-buildver hg38 \
-out results \
-remove \
-protocol refGene,cytoBand \
-operation g,r

# Pull out only exonic or splicing
awk -F'\t' '$6=="exonic" || $6=="splicing" || NR==1' \
results.hg38_multianno.txt > results.exonic_splicing.txt