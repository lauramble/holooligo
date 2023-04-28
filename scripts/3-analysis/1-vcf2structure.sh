VCF="FUT1.vcf"
NAME=${$VCF%.*}
PLINK="~/plink/plink"

$PLINK --vcf $VCF --recode-structure --out $NAME --allow-extra-chr
awk 'NR>2' $NAME.recode.strct_in > $NAME.recode.noheader.strct_in
