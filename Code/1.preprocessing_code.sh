# This script is meant to be a guide for the methodology used and the order. 
# It is not meant to be a standalone code that can work for all workflows. 
# This script also requires bowtie2 indices to be built for the negative controls (swab and blank).
# This script uses conda environmemts to run packages. The following packages could all be installed into the same environment, except for Zebra. 
# Zebra is a standalone set of python scripts. Refer Hakim at al. 2022 for more information.


#This script uses the following packages:
#
#1. Trimgalore (Martin, 2011)
#2. Flash2 (Magoc and Salzberg, 2011)
#3. FASTQC (Babraham Bioinformatics, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#4. MultiQC (Seqera Labs, https://docs.seqera.io/multiqc) 
#5. Bowtie2 (Langmead et al., 2012)
#6. woltka (Zhu et al., 2022)
#7. Zebra (Hakim et al., 2022)

#This script also requires the following dataset to be downloaded
#1. Web Of Life (Zhu et al., 2019)
#2. Skin Microbial Genome Collection (SMGC) (Saheb-Kashaf et al. 2022)

#--------------------------------------------------------------------


#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

conda activate woltka

#making directories

mkdir flash2_output
mkdir fastqc_reports
mkdir merged_flash
mkdir samfiles
mkdir Classify_Species_zebra_25perc

#sample names:

ls -1 *.fastq.gz|rev|cut -d '_' -f 1,2,3 --complement|rev|sort|uniq > sample_names.txt

#loop for trimming, fastqc and flash2 merge

while read -r sample; do
	trim_galore -q 20 -j 4 --paired ./${sample}R1_001.fastq.gz ./${sample}R2_001.fastq.gz -fastqc
	mv *fastqc* ./fastqc_reports
	flash2 ./${sample}R1_001_val_1.fq.gz ./${sample}R2_001_val_2.fq.gz -x 0.10 -M 150 -t 16 -d ./flash2_output --suffix ${Sample}flash.fq.gz| cat > ./merged_flash/${sample}merged.fq.gz
done <sample_names.txt

#QC using MultiQC

multiqc ./fastqc_reports/

# (optional) Computing average read length within each sample. 

while read -r sample; do
        average=$(awk 'NR >= 2 && (NR - 2) % 4 == 0 { char_count[length]++; total_chars += length; } END { avg = total_chars / (((NR - 2) / 4) + 1); print "Average: " avg }' "${sample}merged.fq.gz")
        echo "${sample}\t${average}"   >> ./avg_read_lengths.csv
        echo "${sample} done"
done < sample_names.txt

# Make bowtie2 indices for swab and blank samples. 

bowtie2-build swab_blank_control_30uLbeads.fasta .../swab_blank_control_index


# filtering steps: from each experimental file, filter out swab, blank and host reads. Swab and blank filtered out with global mode. Host filtered out with local. 

while read -r sample; do

        echo -e "\nfiltering swab and blank reads from ${sample} files in paired mode"

        bowtie2 -p 10 -x .../swab_blank_control_index/swab_blank_control_index \
                -1 .../trimmed_files_backup/${sample}_L002_R1_001_val_1.fq.gz \
                -2 .../trimmed_files_backup/${sample}_L002_R2_001_val_2.fq.gz \
                --very-sensitive \
                --un-conc-gz ${sample}_swabandblanks_removed \
                > /dev/null 2> "${sample}_swabblanks_metrics.txt"

        mv ${sample}_swabandblanks_removed.1 .../swab_blank_removed/${sample}_swabandblanks_removed_1.fastq.gz
        mv ${sample}_swabandblanks_removed.2 .../swab_blank_removed/${sample}_swabandblanks_removed_2.fastq.gz

        echo -e "\ndone filtering swab and blank reads for ${sample} files"

        SALGN=$(grep "overall alignment rate" ${sample}_swabblanks_metrics.txt| awk '{print $1}')

        echo -e "\nfiltering host reads from ${sample} files in paired mode"

        bowtie2 -p 8 -x .../databases/GRCh38p14/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index \
                -1 .../swab_blank_removed/${sample}_swabandblanks_removed_1.fastq.gz \
                -2 .../swab_blank_removed/${sample}_swabandblanks_removed_2.fastq.gz \
                --very-sensitive-local \
                --un-conc-gz ${sample}_host_removed \
                > /dev/null 2> "${sample}_host_metrics.txt"

        mv ${sample}_host_removed.1 .../host_removed/${sample}_host_removed_1.fastq.gz
        mv ${sample}_host_removed.2 .../host_removed/${sample}_host_removed_2.fastq.gz

        echo -e "\ndone filtering host reads for ${sample} files"

        HALGN=$(grep "overall alignment rate" ${sample}_host_metrics.txt| awk '{print $1}')

        echo -e "\ndone filtering host reads from all files"

        echo -e "\nmoving on to SMGC alignment rates"

        echo -e "\naligning ${sample} host-removed files to SMGC database"

        bowtie2 -p 8 -x .../databases/SMGC/bacteria_mags_nocontigs/bacteria_mags_nocontigs \
                -1 .../host_removed/${sample}_host_removed_1.fastq.gz \
		-2 .../host_removed/${sample}_host_removed_2.fastq.gz \
                --very-sensitive-local \
                --met-file "${sample}_SMGC_metrics.txt" \
                -S /dev/stdout | gzip > "${sample}_SMGC.sam.gz"

        DALGN=$(grep "overall alignment rate" ${sample}_SMGC_metrics.txt| awk '{print $1}')

        echo -e "${sample}\t${SALGN}\t${HALGN}\t${DALGN}" >> host_alignment_stats.tsv

 	echo -e "\ndone. aligning to WoL"

   	bowtie2 -U ./merged_flash/${sample}merged.fq.gz -x ./*location_of_WoL*/WoLr1 
    		-p 12 
      		--very-sensitive --no-head --no-unal -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" | cut -f1-9 | sed 's/$/\t*\t*/' | gzip > ./WoL_samfiles/${sample}.sam.gz

done <sample_names.txt

echo -e "\ncompleted filtering all samples and then aligning to SMGC. Refer host_alignment_stats.tsv for filtering and alignment details"


# Calculate genomes coverages (onbly for SMGC alignment. WoL was used to identify relic-DNA for beta-diversity only. If using for taxonomy, must do coverage calculations)

python .../zebra_filter/calculate_coverages.py -i ./samfiles -o SMGC_coverages.txt -d ~/storage/databases/SMGC/SMGC_bacteria_metadata.tsv

# Move to R script ("Rmd/0.Zebra_coverages_calculated_aligningtoSMGC") to determine what coverage threshold top pick. This should be the global minima in the genome coverage density plots. 








