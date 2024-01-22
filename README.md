# ITLN1_paper
This repository contains processed datasets and code for generating figures in the Nature Communications paper, _A common polymorphism in the Intelectin-1 gene influences mucus plugging in severe asthma_.

All scripts are in the 'Code' directory and data referenced in the scripts are in the 'Data' directory. 

Here is a summary of where gene expression/protein count data as well as metadata for all samples are located:

**Large bulk RNA-seq datasets**
Large count datasets, must be downloaded from GEO:
1. [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145013) for a Seurat object containing the bronchial brushing single cell data along with cell type annotations for two donors (Z002KB_fresh = Female; Z003TS_fresh = Male) (Seurat object = Z002KB_3TS_integrated_culled.rds).
2. [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152004) for *in vivo* bulk sequence data for GALA II (raw counts = GSE152004_695_raw_counts.txt.gz).

Relevant data from RNA-seq datasets and sample metadata can be found in the Data directory (as excel files):
1. Data/phen.xlsx - for *in vivo* bulk sequence data for GALA II.
2. Data/v2_t2.df.xlsx - for *in vivo* bulk sequence data for SARP.
3. Data/Tracheal_donors_40_metadata.xlsx - for *in vitro* ITLN1 network eigengenes and sample metadata for 40 tracheal donors.

**Human bronchial epithelial culture (HBEC) datasets**
Count datasets from IL-13-stimulated human bronchial epithelial cultures (HBECs) can be found in the Data directory (as text files):
1. Data/secRna_353.txt - Ampliseq RNA-seq counts from 19 donors.
2. Data/secretome_count_matrix.txt - Peptide counts from the aqueous fraction of apical ALI secretions from 14 donors.
3. Data/secretome_count_matrix_newMucus.txt - Peptide counts from the mucus fraction of apical ALI secretions from 9 donors.

Metadata for HBEC donors can be found in the Data directory (as a text file).
1. Data/secretome_metadata.txt

**Human tracheal epithelial culture datasets**
Metadata for tracheal donors used in the mucociliary and ciliary beat freqency analyses are in the Data directory:
1. Data/Tracheal_donor_metadata.xlsx
