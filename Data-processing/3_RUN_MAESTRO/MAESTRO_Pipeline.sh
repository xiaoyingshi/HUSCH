# Note:
# The analysis pipeline starts from a count matrix, which can be MTX, HDF5 or plain text formatted, and then perform QC, clustering and cell-type annotation.
# The analysis pipeline generates a R script which can be modified, clustering and annotated UMAP plots and a Seurat object which can be readable and writable.
# Meta data is optional, not required.

# activate conda environment
source /fs/home/shixiaoying/miniconda3/envs/MAESTRO_normal


# run the following command for detailed usage of MAESTRO
MAESTRO --help
MAESTRO scrna-analysis --help

# the following commands can be used for format conversion
MAESTRO mtx-to-h5           # Convert 10X mtx format matrix to HDF5 format.
MAESTRO mtx-to-count        # Convert 10X mtx format matrix to plain text count table.
MAESTRO count-to-h5         # Convert plain text count table to HDF5 format.
MAESTRO h5-to-count         # Convert HDF5 format to plain text count table.
MAESTRO merge-h5            # Merge 10X HDF5 files.

# Note:
# If there're multiple samples (count matrices) in one dataset, usually we convert the count matrices to HDF5 format firstly. 
# Then we merge these h5 matrices to a big h5 matrix and take it as the input of the analysis pipeline.
# If meta table is provided, users should make sure the barcodes in the meta data are consistent with those in the count matrix.
# Because in the big merged matix, the cell barcodes may have prefix (added through the parameter --cellprefix in MAESTRO merge-h5).


# Examples:
# run analysis pipeline (h5 data)
MAESTRO scrna-analysis --mode normal --format h5 --matrix fileName.h5 \
--meta-file meta.tsv --meta-sep tab --meta-cell 1 \
--count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 \
--directory . --outprefix fileName

# run analysis pipeline (mtx data)
MAESTRO scrna-analysis --mode normal --format mtx --matrix fileName.mtx  \
--meta-file meta.tsv --meta-sep tab --meta-cell 1 \
--feature genes.tsv --gene-column 2 --barcode barcodes.tsv \
--count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 \
--directory . --outprefix fileName

# run analysis pipeline (plain text)
MAESTRO scrna-analysis --mode normal --format plain --matrix fileName.txt --separator tab \
--meta-file meta.tsv --meta-sep tab --meta-cell 1 \
--count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 \
--directory . --outprefix fileName


# Note:
# All the parameters used in example code are set for normal tissue single-cell RNA-seq data.