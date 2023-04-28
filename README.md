# HoloOligo WP2 task 3



## Données

Les données déjà traitées et les fichiers nécessaires au traitement de nouvelles données se trouvent sur le cluster genotoul dans le dossier **`/work/genphyse/genepi/holooligo`**

| Fichier | Description |
|--|--|
| all.bed | Fichier bed définissant des régions d'1Mb autour de gènes d'intérêt : <br> --- 9 gènes FUT <br> --- QTL1, GRP107, MC1R, SOX9, ASIP |
|corrected_Sus_scrofa.Sscrofa11.1.dna.toplevel.fa|Ficher fasta corrigé pour une annotation correcte pour FUT1 (insertion nucléotide **G** à la position **6:59,079,640**)|
|futs.sorted.gtf.gz|Fichier d'annotation gtf corrigé pour une annotation correcte de FUT1 |
|futonly.bed|Fichier bed couvrant les gènes FUT avec 1kb en amont et en aval|
|sus_scrofa_2.vcf|dbSNP modfié pour retirer les espaces|
|gvcf-all|Dossier consitutant la base de donnée GATK (générée avec la commande [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/4404607512475-GenomicsDBImport))|

## Dépendances

Liste des logiciels utilisés :
- **GATK** *(4.2.6.1)*
- **R** *(4.2.2)*
- **bcftools** *(1.11)*
- **Beagle** *(4.0)*
- **perl** *(5.32)*
- **PHASE** *(2.1.1)*
- **VEP** *(109)*

Les packages R nécessaires peuvent être installés avec le script *requirements.R*.

# Utilisation

## Préparation des données

Ces étapes permettent de mettre à jour la [base de données GATK](https://gatk.broadinstitute.org/hc/en-us/articles/4404607512475-GenomicsDBImport) et de filtrer, phaser et annoter les vcfs.

### Étape 0 : Mise à jour de la base de données 
#### Dossier de travail : *0-import_to_db*
#### Entrée : 
- SAMPLE_MAP : voir fichier exemple ou [documentation GATK](https://gatk.broadinstitute.org/hc/en-us/articles/4404607512475-GenomicsDBImport)
- REGION : région cible sous forme de chaîne de caractères (ex: 1:0-10000) ou chemin vers fichier *.bed*
- DB : chemin vers la base de données
#### Procédure :
- Ouvrir le script `import_to_db.sh`et adapter les variables d'entrées aux données
- Lancer le script `import_to_db.sh` sur genotoul (ou copier-coller dans terminal genotoul)
#### Sortie :
- Aucune

### Etape 1 : Extraction vcf

#### Dossier de travail: 1-extract_from_db

#### Entrée :
- REF : génome de référence
- DB : chemin vers la base de données
- REGION : région cible sous forme de chaîne de caractères (ex: 1:0-10000) ou chemin vers fichier *.bed*
- NAME : nom du fichier de sortie
- MINQUAL : seuil de qualité pour filtrer les variants
- DBSNP : chemin vers le fichier dbSNP utilisé

#### Sortie :
- fichiers vcf intermédiaires
- fichier *<NAME\>.vcf* final

#### Procédure :
- Ouvrir le script `extract_from_db.sh`et adapter les variables d'entrées aux données
- Lancer le script `extract_from_db.sh` sur genotoul (ou copier-coller dans terminal genotoul)

### Etape 2 : Phasage et annotation

#### Dossier de travail: 2-phase_annotation

#### Entrée :
- FILE :
- PATH_TO_VEP : chemin vers VEP
- PATH_TO_PHASE : chemin vers PHASE

#### Sortie :
- fichier VCF phasé *.phased.vcf* 
- fichier VCF phasé et annoté *.vep.vcf*
- fichiers de sortie PHASE et VEP

#### Procédure :
- Ouvrir le script `run_phase_annotation.sh`et adapter les variables d'entrées aux données
- Lancer le script `run_phase_annotation.sh` (ou copier-coller dans terminal)

## Analyses et figures
Ces étapes permettent de générer des figures et tables pour faciliter l'analyse des données. Voir section suivante pour une explication détaillée de chaque script

### Etape 3 : Analyse

#### Dossier de travail: 3-analysis

#### Entrée :
- VCF : fichier *.vcf* d'entrée
- POPULATION : fichier *.tsv*

#### Sortie :
- figures *.png*
- tables *.tsv*
- pages *.html*

#### Procédure :
- Ouvrir le script `1-vcf2structure.sh`et adapter les variables d'entrées aux données
- Lancer le script `1-vcf2structure.sh` (ou copier-coller dans terminal)
- Ouvrir le script `2-run_analysis.sh`et adapter les variables d'entrées aux données
- Lancer le script `2-run_analysis.sh` (ou copier-coller dans terminal)

# Scripts R
Les analyses sont regroupées dans plusieurs scripts R. Il est possible de les lancer en ligne de commande (voir script `2-run_analysis.sh`) ou de les ouvrir dans RStudio afin de pouvoir ajuster le script. Pour ce faire, il suffit de commenter les sections indiquées et de dé-commenter les lignes contrôlant les variables (également indiquées dans le script).
## preprocessing.R

Permet de générer les fichiers nécessaires aux scripts suivants ainsi que de donner une vue d'ensemble sur les données.

#### Sortie :
- *variant_table.tsv* : liste des variants annotés
- *reference_allele_freq.tsv*: fréquence de l'allèle de référence pour chaque variant
- *population_overview.png* : distribution des populations
- *population_breakdown.tsv* : tableau donnant le nombre d'individus par population
- *reference_allele_freq_pop.tsv*, *reference_allele_freq_pop.png* : fréquence de l'allèle de référence pour chaque variant par population

## PCA.R
Permet de générer des PCA d'ensemble et par population, ainsi que des liste de SNP corrélées

#### Sortie :
- *ggplot* : figures PCA 2D générées avec ggplot
- *plotly* : figures PCA 2D générées avec plotly
- *plotly3d* : figures PCA 3D générées avec plotly
- *bypop* : figures PCA par population
- *selected_snp.tsv* : tableau des SNP corrélés avec les 4 premières composantes principales

## structure.R
Permet de générer des figures de type Structure

#### Sortie :
- *cross_entropy.png* : cross-entropy graph to choose the best number of clusters

## haplotypes.R

Permet d'étudier les haplotypes

#### Sortie :
- *haplotype_number_bypop.tsv*, *haplotype_number.tsv* : quantités d'haplotypes (ensemble et par population)
- *allele_frequency_majority.png*, *allele_frequency_all.png* : fréquence allélique de l'allèle de référence dans les haplotypes
-*allele_frequency_majority.png*, *haplotype_distribution_all.png* : distribution des haplotypes par population

## combine.R
Permet de visualiser les combinaisons entre haplotypes
#### Sortie :
Aucune

