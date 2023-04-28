# RUN ON GENOTOUL

module load bioinfo/Beagle_v4.0
module load bioinfo/gatk-4.2.6.1
module load bioinfo/bcftools-1.11

REF="/work/genphyse/genepi/genomes/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
DB="/work/genphyse/holooligo/gvcf-all"
REGION="6:54076431-54081514" # Also accepts .bed files
NAME="FUT1"
MINQUAL="20.0"
DBSNP="work/genphyse/genepi/holooligo/reference/sus_scrofa_2.vcf.gz"

srun --mem 16g gatk --java-options "-Xmx16g" GenotypeGVCFs \
	-R  $REF \
	-V gendb://$DB \
	-L $REGION \
	-O $NAME.vcf.gz \
	--standard-min-confidence-threshold-for-calling $MINQUAL \
	--allow-old-rms-mapping-quality-annotation-data &> genotypegvcfs.log

srun bcftools view -M2 -m2 $NAME.vcf.gz -o $NAME.filtered.vcf.gz -O z # filter biallelic
tabix $NAME.filtered.vcf.gz

srun bcftools annotate -a $DBSNP -c ID -o $NAME.annotated.filtered.vcf.gz -O z $NAME.filtered.vcf.gz # annotate ID (rs00000)

srun --mem=32G java -Xmx32g -jar $BEAGLE gt=$NAME.annotated.filtered.vcf.gz out=$NAME.annotated.filtered.beagle idb=true &> beagle.log

srun bgzip -d $NAME.annotated.filtered.beagle.vcf.gz

mv $NAME.annotated.filtered.beagle.vcf $NAME.vcf
