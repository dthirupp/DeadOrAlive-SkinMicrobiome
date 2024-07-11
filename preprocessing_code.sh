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

#--------------------------------------------------------------------


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

#Aligning to Web Of Life

while read -r sample; do
	bowtie2 -U ./merged_flash/${sample}merged.fq.gz -x ./*location_of_WoL*/WoLr1 -p 12 --very-sensitive --no-head --no-unal -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" | cut -f1-9 | sed 's/$/\t*\t*/' | gzip > ./samfiles/${sample}.sam.gz
done < sample_names.txt

#Calculate genomes coverages

source ~/.bashrc
conda activate woltka
python ~/zebra_filter/calculate_coverages.py -i ./samfiles -o ./coverages.txt -d ~/zebra_filter/databases/WoL/metadata.tsv


#Create of a subset of the WoL database that contains only genomes with >25% genome coverage:

python ~/zebra_filter/filter_sam.py -i ~/coverages.txt -s ~/*location_of_WoL*/WoLr1 -c 0.25 -o ./Zebra_25perc_samfiles

#Taxonomic Classification at the Genus level

woltka classify --input ./Zebra_25perc_samfiles --map ~/*location_of_WoL*/taxonomy/taxid.map --nodes ~/*location_of_WoL*/taxonomy/nodes.dmp --names ~/*location_of_WoL*/taxonomy/names.dmp --output ./Classify_Species_zebra_25perc/Skin_LiveDead_Zebra_25perc_Woltka_Classify_Species.biom --rank species --name-as-id --outmap ./Classify_Species_zebra_25perc

#Functional annotation with read-to-genus maps for stratification

woltka classify --input ./Zebra_25perc_samfiles --coords ~/*location_of_WoL*/annotation/coords.txt.xz -m ~/*location_of_WoL*/annotation/uniref.map.xz -m ~/*location_of_WoL*/annotation/go/process.map.xz --map-as-rank --rank process -n ~/*location_of_WoL*/annotation/uniref.names.xz -r uniref --stratify ./Classify_Species_zebra_25perc --output ./Classify_Species_zebra_25perc/Skin_LiveDead_Zebra_25perc_Classify_Species_Uniref.tsv

#making OGU tables for core-metrics

pip install upgrade numpy ##upgrade numpy if needed

woltka classify -i ./Zebra_25perc_samfiles -o ogu_25percent_table.biom



conda deactivate








