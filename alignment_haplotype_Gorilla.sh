#!/bin/bash 

#SBATCH -J calling_variants_and_multihetsep
#SBATCH -A SCALLY-SL3-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60000
#SBATCH --get-user-env
#SBATCH --array=1-23
#SBATCH --mail-type=NONE
#SBATCH --time=01:00:00
#SBATCH --error=/home/jg2162/rds/hpc-work/log_files_meta/log_file_species_new_mask/Gorilla/calling_variants%j.err
#SBATCH --output=/home/jg2162/rds/hpc-work/log_files_meta/log_file_species_new_mask/Gorilla/calling_variants%j.out

CHROMS_1=( "chr1_pat_hsa1" "chr2_pat_hsa3" "chr3_pat_hsa4" "chr4_pat_hsa17x5" "chr5_pat_hsa6" "chr6_pat_hsa7" "chr7_pat_hsa8" "chr8_pat_hsa10" "chr9_pat_hsa11" "chr10_pat_hsa12" "chr11_pat_hsa2b" "chr12_pat_hsa2a" "chr13_pat_hsa9" "chr14_pat_hsa13" "chr15_pat_hsa14" "chr16_pat_hsa15" "chr17_pat_hsa18" "chr18_pat_hsa16" "chr19_pat_hsa5x17" "chr20_pat_hsa19" "chr21_pat_hsa20" "chr22_pat_hsa21" "chr23_pat_hsa22" )
CHROMS_2=( "chr1_mat_hsa1" "chr2_mat_hsa3" "chr3_mat_hsa4" "chr4_mat_hsa17x5" "chr5_mat_hsa6" "chr6_mat_hsa7" "chr7_mat_hsa8" "chr8_mat_hsa10" "chr9_mat_hsa11" "chr10_mat_hsa12" "chr11_mat_hsa2b" "chr12_mat_hsa2a" "chr13_mat_hsa9" "chr14_mat_hsa13" "chr15_mat_hsa14" "chr16_mat_hsa15" "chr17_mat_hsa18" "chr18_mat_hsa16" "chr19_mat_hsa5x17" "chr20_mat_hsa19" "chr21_mat_hsa20" "chr22_mat_hsa21" "chr23_mat_hsa22" )

TEMP="/home/jg2162/rds/hpc-work/temp/Gorilla"
OUTPUT="/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Data_species_specific/Gorilla/Output_vcf_file_haplotype"
OUTPUT_MULTIHET="/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Data_species_specific/Gorilla/Output_multihet_from_haplotype"

export PATH=$PATH:/home/jg2162/Cam_simulation_Data/tools/minimap2

CHROM1=${CHROMS_1[$SLURM_ARRAY_TASK_ID-1]}
CHROM2=${CHROMS_2[$SLURM_ARRAY_TASK_ID-1]}

HAP1="/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Diploid_reference_species/Gorilla/mGorGor1.pat.cur.20231122.fasta.gz"
HAP2="/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Diploid_reference_species/Gorilla/mGorGor1.mat.cur.20231122.fasta.gz"

OUTPUT_BED=/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Data_species_specific/Gorilla/Output_BED_file/${CHROM1}
mkdir -p ${OUTPUT_BED}

OUTPUT_INDEL_BED=/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Data_species_specific/Gorilla/Output_indel_BED_file/${CHROM1}
mkdir -p ${OUTPUT_INDEL_BED}

OUTPUT_BAM_CALLER_BED=/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Data_species_specific/Gorilla/Output_bam_caller_BED_file/${CHROM1}
mkdir -p ${OUTPUT_BAM_CALLER_BED}

OUTPUT_MULTIHET_FILTERED_95=/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Data_species_specific/Gorilla/Output_multihet_filtered_95
mkdir -p ${OUTPUT_MULTIHET_FILTERED_95}

OUTPUT_MHS_BOOTSTRAPPED=/home/jg2162/rds/hpc-work/Cam_Evan_eichler_new_data/Data_species_specific/Gorilla/Output_Multihet_Bootstrapped
mkdir -p ${OUTPUT_MHS_BOOTSTRAPPED}

GENERATE_MULTIHETSEP="/home/jg2162/Cam_simulation_Data/tools/msmc-tools/generate_multihetsep.py"
MASK_FILE_GENERATOR="/home/jg2162/Cam_simulation_Data/script/mask_file.py"
BOOTSTRAPPING_GENERATOR="/home/jg2162/Cam_simulation_Data/tools/cobraa/reproducibility/block_bootstrap_mhs.py"
BAM_CALLER="/home/jg2162/Cam_simulation_Data/tools/msmc-tools/bamCaller.py"
INDEL_MASKER="/home/jg2162/Cam_simulation_Data/script/indel_mask.py"

##Extracting the relevant FASTA sequences
# echo "Extracting the chromosome FASTA sequences .... "

samtools faidx $HAP1 $CHROM1 > ${TEMP}/${CHROM1}.fasta
samtools faidx $HAP2 $CHROM2 > ${TEMP}/${CHROM2}.fasta

# echo "Extracting  the chromosome FASTA sequences. Done"

samtools faidx ${TEMP}/${CHROM1}.fasta
samtools faidx ${TEMP}/${CHROM2}.fasta
# ##Mapping using the minimap

# echo "Aligning the chromosome using minimap ...."
minimap2 -ax asm5 ${TEMP}/${CHROM1}.fasta ${TEMP}/${CHROM2}.fasta > ${TEMP}/${CHROM1}.sam
samtools sort ${TEMP}/${CHROM1}.sam -o ${TEMP}/${CHROM1}_sorted.sam
samtools view -S -b ${TEMP}/${CHROM1}_sorted.sam -o ${TEMP}/${CHROM1}.bam


echo "Aligning the chromosome using minimap. Done"

##Indexing the BAM file

echo "Indexing the BAM file ...."
samtools index ${TEMP}/${CHROM1}.bam

echo "Indexing the BAM file. Done"

##mpileup and call the variants

echo "Piling up the variants ...."
bcftools mpileup -f ${TEMP}/${CHROM1}.fasta ${TEMP}/${CHROM1}.bam -Oz -o ${TEMP}/${CHROM1}.vcf.gz
python ${INDEL_MASKER} --input_gen_vcf ${TEMP}/${CHROM1}.vcf.gz --Chromosome_name ${CHROM2} --output_BED ${OUTPUT_INDEL_BED}/${CHROM1}_indel_bed.txt
bcftools call --ploidy 2 -c -V indels ${TEMP}/${CHROM1}.vcf.gz | ${BAM_CALLER} 1 ${OUTPUT_BAM_CALLER_BED}/${CHROM2}_bed.txt.gz | gzip -c > ${TEMP}/${CHROM1}_only_snps.vcf.gz

rm  ${TEMP}/${CHROM1}.vcf.gz

zcat ${OUTPUT_BAM_CALLER_BED}/${CHROM2}_bed.txt.gz | \
awk -v chrom2="$CHROM2" 'BEGIN{OFS="\t"} { $1=chrom2; print $0 }' | \
gzip -c > ${OUTPUT_BAM_CALLER_BED}/${CHROM2}_modified_bed.txt.gz

tabix -p vcf ${TEMP}/${CHROM1}_only_snps.vcf.gz

zcat ${TEMP}/${CHROM1}_only_snps.vcf.gz | sed  's/1\/1/1\/0/g' | bgzip -c > ${TEMP}/${CHROM1}_restructured.vcf.gz

echo "Piling up the variants. Done"

bcftools view -i 'INFO/DP=1' -O z ${TEMP}/${CHROM1}_restructured.vcf.gz -o ${TEMP}/${CHROM1}_restructured_DP_1.vcf.gz

echo "Removing Indels and missing data from vcf file ...."
vcftools --gzvcf ${TEMP}/${CHROM1}_restructured_DP_1.vcf.gz \
         --remove-indels \
         --recode \
         --out ${OUTPUT}/${CHROM1}

echo "Removing Indels and missing data from vcf file. Done"

echo "Bgziping the vcf file ...."
bgzip ${OUTPUT}/${CHROM1}.recode.vcf

echo "Bgziping the vcf file. Done"

echo "Removing the Temporary Directories..."

python $GENERATE_MULTIHETSEP ${OUTPUT}/${CHROM1}.recode.vcf.gz > ${OUTPUT_MULTIHET}/${CHROM1}_multihet.txt

echo "Removing the Temporary Directories, Done. "

echo "Generating the mask for multihet file"
python $MASK_FILE_GENERATOR --input_file ${OUTPUT_MULTIHET}/${CHROM1}_multihet.txt --Chromosome_name $CHROM2 --output_BED $OUTPUT_BED

echo "Mask generated"

echo "Generating the multihetsep file" 

python $GENERATE_MULTIHETSEP ${OUTPUT}/${CHROM1}.recode.vcf.gz --negative_mask="${OUTPUT_BED}/${CHROM2}_95.txt" --mask="${OUTPUT_BAM_CALLER_BED}/${CHROM2}_modified_bed.txt.gz" --negative_mask="${OUTPUT_INDEL_BED}/${CHROM1}_indel_bed.txt" > ${OUTPUT_MULTIHET_FILTERED_95}/${CHROM1}_multihet.txt

echo "The main multihet file generated"

echo "Bootstrapping the multihet file"

for i in $(seq 0 29); do 
   OUTPUT_MHS_BOOTSTRAPPED_I=${OUTPUT_MHS_BOOTSTRAPPED}/${i}
   mkdir -p ${OUTPUT_MHS_BOOTSTRAPPED_I}
   python $BOOTSTRAPPING_GENERATOR --inmhs ${OUTPUT_MULTIHET_FILTERED_95}/${CHROM1}_multihet.txt --windowsize 5e+06 --outmhs ${OUTPUT_MHS_BOOTSTRAPPED_I}/${CHROM1}_multihet.txt &
done
wait
