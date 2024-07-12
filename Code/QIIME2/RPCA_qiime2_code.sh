#Using qiime2 (qiime2-2022.11) RPCA feature via DEICODE plugin

conda activate qiime2-2022.11

#1. Computing RPCA on all samples to observe the effect of PMA treatment on the samples as whole (Fig S1)

#convert tsv file to biom

	biom convert -i Final_rpca_all.tsv -o taxa_all_table.biom --to-hdf5
	
#import file as qza qiime2 artifact 
	
	qiime tools import --type FeatureTable[Frequency] --input-path taxa_all_table.biom --output-path taxa_all_table.qza

#run gemelli auto-rpca to compute robust aithcisson distance

	qiime gemelli auto-rpca --i-table taxa_all_table.qza --p-min-feature-count 10 --p-min-sample-count 500 --o-biplot taxa_all_ordination.qza --o-distance-matrix taxa_all_distance.qza
	
#export the qza artifact back to tsv to use for plotting with ggplot2 in R.

	qiime tools export --input-path taxa_all_distance.qza --output-path taxa_all_distance
	qiime tools export --input-path taxa_all_ordination.qza --output-path taxa_all_ordination
	
	mv /taxa_all_ordination/ordination.txt /taxa_all_ordination.txt

#cleanup

	rm -r /taxa_all_ordination/

#compute biplot

	qiime emperor biplot --i-biplot ./Manuscript_RA_stuff/taxa_all_ordination.qza --m-sample-metadata-file metadata.tsv --o-visualization ./Manuscript_RA_stuff/taxa_all_biplot.qzv --p-number-of-features 2

#calculate significance of treatment on sample clustering using PERMANOVA

	qiime diversity beta-group-significance --i-distance-matrix ./Manuscript_RA_stuff/taxa_all_distance.qza --m-metadata-file metadata.tsv --m-metadata-column PMA_treated --p-method permanova --o-visualization ./Manuscript_RA_stuff/taxa_all_PMA.qzv


	mv /taxa_all_distance/*.tsv ./taxa_all_distance.tsv
	rm -r ./Manuscript_RA_stuff/taxa_all_distance

#viewing of the beta-group significance plot
	qiime tools view ./Manuscript_RA_stuff/taxa_all_PMA.qzv


#2. Computing RPCA on only Raw samples and evaluate the significance of Skin Site Type and Individual on clustering of samples (Fig 2)

#convert tsv file to biom

	biom convert -i Final_rpca_raw.tsv -o taxa_nopma_table.biom --to-hdf5

#import file as qza qiime2 artifact 

	qiime tools import --type FeatureTable[Frequency] --input-path taxa_nopma_table.biom --output-path taxa_nopma_table.qza

#run gemelli auto-rpca to compute robust aithcisson distance

	qiime gemelli auto-rpca --i-table taxa_nopma_table.qza --p-min-feature-count 5000 --p-min-sample-count 500 --o-biplot taxa_no_pma_ordination.qza --o-distance-matrix taxa_nopma_distance.qza

#export the qza artifact back to tsv to use for plotting with ggplot2 in R.

	qiime tools export --input-path taxa_nopma_distance.qza --output-path taxa_nopma_distance
	qiime tools export --input-path taxa_no_pma_ordination.qza --output-path taxa_no_pma_ordination

	mv taxa_no_pma_ordination/ordination.txt taxa_no_pma_ordination.txt

#cleanup

	rm -r /taxa_no_pma_ordination/

#compute biplot

	qiime emperor biplot --i-biplot taxa_no_pma_ordination.qza --m-sample-metadata-file metadata.tsv --o-visualization taxa_nopma_biplot.qzv --p-number-of-features 2

#calculate significance of Skin site type on sample clustering using PERMANOVA

	qiime diversity beta-group-significance --i-distance-matrix taxa_nopma_distance.qza --m-metadata-file metadata.tsv --m-metadata-column Body_Site_Type --p-method permanova --o-visualization taxa_nopma_BSType.qzv

#calculate significance of Individual on sample clustering using PERMANOVA

	qiime diversity beta-group-significance --i-distance-matrix taxa_nopma_distance.qza --m-metadata-file metadata.tsv --m-metadata-column Volunteer_ID --p-method permanova --o-visualization taxa_nopma_Indiv.qzv

	mv /taxa_nopma_distance/*.tsv taxa_nopma_distance.tsv

#cleanup

	rm -r taxa_nopma_distance


#viewing of the beta-group significance plots

	qiime tools view taxa_nopma_BSType.qzv
	qiime tools view taxa_nopma_Indiv.qzv


#3. Computing RPCA on only PMA-treated samples and evaluate the significance of Skin Site Type and Individual on clustering of samples (Fig 2)

#convert tsv file to biom

        biom convert -i Final_rpca_pma.tsv -o taxa_pma_table.biom --to-hdf5

#import file as qza qiime2 artifact 

        qiime tools import --type FeatureTable[Frequency] --input-path taxa_pma_table.biom --output-path taxa_pma_table.qza

#run gemelli auto-rpca to compute robust aithcisson distance

        qiime gemelli auto-rpca --i-table taxa_pma_table.qza --p-min-feature-count 5000 --p-min-sample-count 500 --o-biplot taxa_pma_ordination.qza --o-distance-matrix taxa_pma_distance.qza

#export the qza artifact back to tsv to use for plotting with ggplot2 in R.

        qiime tools export --input-path taxa_pma_distance.qza --output-path taxa_pma_distance
        qiime tools export --input-path taxa_pma_ordination.qza --output-path taxa_pma_ordination

        mv taxa_pma_ordination/ordination.txt taxa_pma_ordination.txt

#cleanup

        rm -r /taxa_pma_ordination/

#compute biplot

        qiime emperor biplot --i-biplot taxa_pma_ordination.qza --m-sample-metadata-file metadata.tsv --o-visualization taxa_pma_biplot.qzv --p-number-of-features 2

#calculate significance of Skin site type on sample clustering using PERMANOVA

        qiime diversity beta-group-significance --i-distance-matrix taxa_pma_distance.qza --m-metadata-file metadata.tsv --m-metadata-column Body_Site_Type --p-method permanova --o-visualization taxa_pma_BSType.qzv

#calculate significance of Individual on sample clustering using PERMANOVA

        qiime diversity beta-group-significance --i-distance-matrix taxa_pma_distance.qza --m-metadata-file metadata.tsv --m-metadata-column Volunteer_ID --p-method permanova --o-visualization taxa_pma_Indiv.qzv

        mv /taxa_pma_distance/*.tsv taxa_pma_distance.tsv

#cleanup

        rm -r taxa_pma_distance


#viewing of the beta-group significance plots

        qiime tools view taxa_pma_BSType.qzv
        qiime tools view taxa_pma_Indiv.qzv

conda deactivate

