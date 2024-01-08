# ITLN1_paper
This repository contains code for generating figures in the Nature Communications paper, _A common polymorphism in the Intelectin-1 gene influences mucus plugging in severe asthma_.

All scripts are in the 'Code' directory and data referenced in the scripts are in the 'Data' directory. Large count datasets, however must be downloaded from GEO ([here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145013) for bronchial brushing single cell data, [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152004) for *in vivo* bulk sequence data for GALA II, [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152004) for *in vivo* bulk sequence data for SARP, and [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152004) for *in vitro* bulk sequence data for tracheal cultures. For the single cell datasets, subsequent analyses can be run by first loading saved R files containing the relevant datasets (Seurat objects), which can also be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145013).

Count datasets from IL-13-stimulated bronchial ALI cultures can be found in the 'Data' directory:
1. Ampliseq RN-seq counts from 19 donors (secRna_353.txt).
2. Peptide counts from the aqueous fraction of apical ALI secretions from 14 donors (secretome_count_matrix.txt).
3. Peptide counts from the mucus fraction of apical ALI secretions from 9 donors (secretome_count_matrix_newMucus.txt).
