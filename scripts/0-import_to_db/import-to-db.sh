module load bioinfo/gatk-4.2.6.1

SAMPLE_MAP="update-sample-map.txt"
DB="/work/genphyse/holooligo/gvcf-all"
REGION="/work/genphyse/genepi/holooligo/all.bed"

srun --mem=16G gatk --java-options "-Xmx16g -Xms12g" GenomicsDBImport \
      --sample-name-map $SAMPLE_MAP \
      --genomicsdb-update-workspace \
      --path $DB \
      -L $REGION \
      &> log-update.txt &
