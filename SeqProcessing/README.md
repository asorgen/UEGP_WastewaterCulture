
# UEGP Wastewater Culturing Analysis

## Transparency and Reproducibility

Sequences were processed using QIIME2 version: 2019.1.0 on a general use Slurm partition, Orion.

## Steps for sequence processing

**1.) Generate fastq manifest file and secure copy to research cluster**                

See manifest file `fastq_manifest.csv`



**2.) Import manifest into QIIME2**               

_Run time ~21 min using 16 cores from 1 node and 125gb memory_

`qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path fastq_manifest.csv --input-format PairedEndFastqManifestPhred33 --output-path raw-seqs.qza`



**3.) Summarize and view .qza file**               

_Run time ~59 min using 16 cores from 1 node and 125gb memory_

`qiime demux summarize --i-data raw-seqs.qza --o-visualization raw-seqs.qzv`

Inspect raw-seqs.qzv using https://view.qiime2.org/ to determine subsequent quality filtering steps.


**4.) Denoise and filter sequences using DADA2**               

_Run time ~186 hrs using 16 cores from 1 node and 125gb memory_

`qiime dada2 denoise-paired --i-demultiplexed-seqs raw-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 115 --p-trunc-len-r 115 --o-table ASVtable.qza --o-representative-sequences ASVseqs.qza --o-denoising-stats denoising-stats.qza --verbose --p-n-threads 16`

Visualize sequences and tables

`qiime feature-table summarize --i-table ASVtable.qza --o-visualization ASVtable.qzv`

`qiime feature-table tabulate-seqs --i-data ASVseqs.qza --o-visualization ASVseqs.qzv`



**5.) Remove ASVs with frequency <10**               

`qiime feature-table filter-features --i-table ASVtable.qza --p-min-frequency 10 --o-filtered-table filtered-ASVtable.qza`

Visualize table

`qiime feature-table summarize --i-table filtered-ASVtable.qza --o-visualization filtered-ASVtable.qzv`

Filter low frequency ASVs from sequences

`qiime feature-table filter-seqs --i-data ASVseqs.qza --i-table filtered-ASVtable.qza --o-filtered-data filtered-ASVseqs.qza`



**6.) Import Silva database into QIIME2 artifacts**

`qiime tools import --type 'FeatureData[Sequence]' --input-path /REFdatabase/silva_128_99_otus_16S.fasta --output-path /REFdatabase/silva128-otus-99.qza`

`qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path /REFdatabase/silva_128_99_16S_consensus_taxonomy.txt --output-path /REFdatabase/silva128-consensus-taxonomy-99.qza`



**7.) Open reference OTU clustering at 100% sequence identity** 

`qiime vsearch cluster-features-open-reference --i-table filtered-ASVtable.qza --i-sequences filtered-ASVseqs.qza --i-reference-sequences /REFdatabase/silva128-otus-99.qza --p-perc-identity 1.0 --o-clustered-table OR-OTUtable.qza --o-clustered-sequences OR-OTUseqs.qza --o-new-reference-sequences OR-new-ref-OTUseqs.qza --p-threads 16 --verbose`	

`qiime feature-table summarize --i-table OR-OTUtable.qza --o-visualization OR-OTUtable.qzv`



**8.) Identify chimeras**

`qiime vsearch uchime-denovo --i-table OR-OTUtable.qza --i-sequences OR-OTUseqs.qza --output-dir /OR-chimeraID`



**9.) Filter out chimeras from OTU table**

`qiime feature-table filter-features --i-table OR-OTUtable.qza --m-metadata-file /OR-chimeraID/nonchimeras.qza --o-filtered-table OR-OTUtable-nochim.qza`

`qiime feature-table summarize --i-table OR-OTUtable-nochim.qza --o-visualization OR-OTUtable-nochim.qzv`



**10.) Filter out chimeras from sequences**

`qiime feature-table filter-seqs --i-data OR-OTUseqs.qza --m-metadata-file /OR-chimeraID/nonchimeras.qza --o-filtered-data OR-OTUseqs-nochim.qza`



**11.) Remove singletons from OTU table**

`qiime feature-table filter-features --i-table OR-OTUtable-nochim.qza --p-min-frequency 2 --o-filtered-table OR-OTUtable-nochim-nosingle.qza`

`qiime feature-table summarize --i-table OR-OTUtable-nochim-nosingle.qza --o-visualization OR-OTUtable-nochim-nosingle.qzv`

Inspect OR-OTUtable-nochim-nosingle.qzv using https://view.qiime2.org/ to determine total sequence frequency.

*Table summary*

Number of samples: 309

Number of features: 136,765

Total frequency: 241,320,131


**12.) Filter out singletons from sequences**

`qiime feature-table filter-seqs --i-data OR-OTUseqs-nochim.qza --i-table OR-OTUtable-nochim-nosingle.qza --o-filtered-data OR-OTUseqs-nochim-nosingle.qza`



**13.) Remove OTUs <0.01% across all samples**

0.01% of Total frequency = 24,132

`qiime feature-table filter-features --i-table OR-OTUtable-nochim-nosingle.qza --p-min-frequency 24132 --o-filtered-table OR-OTUtable-nochim-nosingle-abun.qza`

`qiime feature-table summarize --i-table OR-OTUtable-nochim-nosingle-abun.qza --o-visualization OR-OTUtable-nochim-nosingle-abun.qzv`

Inspect OR-OTUtable-nochim-nosingle-abun.qzv using https://view.qiime2.org/

*Table summary*

Number of samples: 309

Total OTUs: 1,156

Total OTU counts: 204,547,544

Further inspect OR-OTUtable-nochim-nosingle-abun.qzv using https://view.qiime2.org/ using 'Interactive Sample Detail' tab to determine minimum sample sequence count.



**14.) Filter out super low reads**

`qiime feature-table filter-samples --i-table OR-OTUtable-nochim-nosingle-abun.qza --p-min-frequency 10000 --o-filtered-table final-filtered-OR-OTUtable.qza`

`qiime feature-table summarize --i-table final-filtered-OR-OTUtable.qza --o-visualization final-filtered-OR-OTUtable.qzv`

Inspect final-filtered-OR-OTUtable.qzv using https://view.qiime2.org/

*Table summary*

Number of samples: 303

Total OTUs: 1,156

Total OTU counts: 204,546,845



**15.) Filter sequences that have been removed from tables**

`qiime feature-table filter-seqs --i-data OR-OTUseqs-nochim-nosingle.qza --i-table final-filtered-OR-OTUtable.qza --o-filtered-data final-filtered-OR-OTUseqs.qza`

`qiime feature-table tabulate-seqs --i-data final-filtered-OR-OTUseqs.qza --o-visualization final-filtered-OR-OTUseqs.qzv`



**16.) Generate tree for phylogenetic diversity analyses**

`qiime phylogeny align-to-tree-mafft-fasttree --i-sequences final-filtered-OR-OTUseqs.qza --o-alignment mafft-aligned-OR-seqs.qza --o-masked-alignment masked-mafft-aligned-OR-seqs.qza --o-tree OR-unrooted-tree.qza --o-rooted-tree OR-rooted-tree.qza`

 

**17.) Calculate alpha and beta diversity based on phylogeny**

`qiime diversity core-metrics-phylogenetic --i-phylogeny OR-rooted-tree.qza --i-table final-filtered-OR-OTUtable.qza --p-sampling-depth 11500 --m-metadata-file metadata.tsv --output-dir diversity-metrics`



**18.) Calculate Shannon alpha diversity and observed OTUs**

`qiime diversity alpha --i-table final-filtered-OR-OTUtable.qza --p-metric 'shannon' --o-alpha-diversity shannon-diversity.qza`

`qiime metadata tabulate --m-input-file shannon-diversity.qza --o-visualization shannon-diversity.qzv`

`qiime diversity alpha --i-table final-filtered-OR-OTUtable.qza --p-metric 'observed_otus' --o-alpha-diversity observed-otus.qza`

`qiime metadata tabulate --m-input-file observed-otus.qza --o-visualization observed-otus.qzv`



**19.) Assign taxonomy to OTUs using Silva 128**

`qiime feature-classifier classify-consensus-vsearch --i-query final-filtered-OR-OTUseqs.qza --i-reference-reads /REFdatabase/silva128-otus-99.qza --i-reference-taxonomy /REFdatabase/silva128-consensus-taxonomy-99.qza --o-classification OR-taxonomy.qza`

`qiime metadata tabulate --m-input-file OR-taxonomy.qza --o-visualization OR-taxonomy.qzv`



**20.) Remove all mitochondria/chloroplast identification from taxonomy**

`qiime taxa filter-table --i-table final-filtered-OR-OTUtable.qza --i-taxonomy OR-taxonomy.qza --p-exclude mitochondria,chloroplast --o-filtered-table OR-OTUtable.qza`

`qiime feature-table summarize --i-table OR-OTUtable.qza --o-visualization OR-OTUtable.qzv`



**21.) Make taxonomy barplots**

`qiime taxa barplot --i-table OR-OTUtable.qza --m-metadata-file metadata.tsv --i-taxonomy OR-taxonomy.qza --o-visualization OR-taxonomy-barplot.qzv`



**22.) Export the OTU IDs and sequences**

`qiime tools export --input-path OR-OTUtable.qza --output-path OR-OTUtable-export`

`biom convert --to-tsv -i OR-OTUtable-export/feature-table.biom -o OR-OTUtable-export/OR-OTUfeature-table.tsv`



**23.) Get OTU taxonomy**

Open up OR-taxonomy.qzv in QIIME2 View

Download the metadata.tsv file

`mv ~/Downloads/metadata.tsv UEGP_WastewaterCulture/data/OTU-taxonomy.tsv`



**24.) Download metadata tables**

Open up shannon_vector.qzv in QIIME2 View

`mv ~/Downloads/metadata.tsv UEGP_WastewaterCulture/data/shannon.tsv`

Open up observed_otus_vector.qzv in QIIME2 View

`mv ~/Downloads/metadata.tsv UEGP_WastewaterCulture/data/observed_otus.tsv`

Open up final-filtered-OR-OTUtable.qzv in QIIME2 View

Download Frequency per sample csv 

`mv ~/Downloads/sample-frequency-detail.csv UEGP_WastewaterCulture/data/SEQfreqperSample.csv`



**25.) Download taxonomy tables**

Open up taxonomy-barplot.qzv in QIIME2 View

Download the csv files

`mv ~/Downloads/level-2.csv UEGP_WastewaterCulture/data/taxonomy-barplot`

`mv ~/Downloads/level-3.csv UEGP_WastewaterCulture/data/taxonomy-barplot`

`mv ~/Downloads/level-4.csv UEGP_WastewaterCulture/data/taxonomy-barplot`

`mv ~/Downloads/level-5.csv UEGP_WastewaterCulture/data/taxonomy-barplot`

`mv ~/Downloads/level-6.csv UEGP_WastewaterCulture/data/taxonomy-barplot`

`mv ~/Downloads/level-7.csv UEGP_WastewaterCulture/data/taxonomy-barplot`

