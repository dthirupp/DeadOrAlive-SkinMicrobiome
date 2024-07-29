while read -r sample; do
	average=$(awk 'NR >= 2 && (NR - 2) % 4 == 0 { char_count[length]++; total_chars += length; } END { avg = total_chars / (((NR - 2) / 4) + 1); print "Average: " avg }' "${sample}_merged.fq.gz")
	echo "${sample}\t${average}"   >> ./avg_read_lengths.csv
	echo "${sample} done"
done < /sample_names_final.txt
