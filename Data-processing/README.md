Reference update pipeline
================
### Workflow

![Workflow](https://raw.githubusercontent.com/SELINA-team/SELINA-reference_construction/master/docs/workflow.png)

### Download data

1.  Download datasets from databases. Using python crawler for GEO, and
    manully download for other databases such as EBI, HCA, and Broad.

### Format data

2.  Formatted raw data to **plain/mtx/h5** and check meta information
    from crawler manually.  

<!-- end list -->

  - 2.1 Split into normal or tumor dataset.  
  - 2.2 Collect meta information(Sample, Patient, Platform, CellType (if has
    original annotation or curated.assign with corresponding markers
    mentioned in paper)  
  - 2.3 Some common formatting examples that may be helpful to you are
    listed in `2_Format_Data/tips.txt`.

### Generate rds file

3.  Run MAESTRO (Please check and fill in ‘Tissue’ and ‘Stage’ columns
    in excel-table)  
    **code: **`3_RUN_MAESTRO/MAESTRO_Pipeline.sh`

### Quality Control

4.  Quality Control  
- 4.1 File preparation: meta file  

<!-- end list -->

  - Format

| rds\_path | batchName |
| --------- | --------- |
| r1        | b1        |
| r2        | b2        |

- 4.2 Check and remove batch effect.

``` bash
Rscript 4_QC/qc.R -i meta_file_path -o output_path -t 8 
```

### Annotation

5.  Annotate single cell data based on markers collected.  
- 5.1 File preparation: Marker file  

<!-- end list -->

  - Format

| Celltype | Marker.gene |
| -------- | ----------- |
| c1       | g1,g2,g3…   |
| c2       | …           |
  
- 5.2 Generate signature list  
- 5.3 Annotate with function  
- 5.4 Use featureplot to check marker expression and curation. If some cluster’s annotation are weird, replace it with the right one.  
- 5.5 Unify into different level
