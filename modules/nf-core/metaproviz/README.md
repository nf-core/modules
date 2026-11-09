# MetaProViz modules

These modules are the building blocks of the upcoming MetaProViz Nextflow pipeline for the pre-processing of metabolomics data. They wrap functions of the R package [MetaProViz](https://saezlab.github.io/MetaProViz/).

## The modules

| Module                        | What it does                                                                                                                                                                         |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `metaproviz/internalstandard` | Quality control with the internal standard (IS) spike-ins. The CV of every IS across the samples shows the technical noise of the measurement, and the IS are removed from the data. |

Further modules will follow, for example `poolestimation`, `processingfilter`, `processingimpute`, `processingnormalize` and `processingoutliers`.

## Input and output format

The input is either a SummarizedExperiment (SE) saved as RDS (`se_rds`) or three TSV or CSV tables (`data_matrix`, `feature_matrix` and `sample_matrix`), and only one of the two shapes is given. The assay of the SE corresponds to the `data_matrix`, the `rowData` to the `feature_matrix` and the `colData` to the `sample_matrix`. In the tables, the first column holds the feature or sample IDs.

The value input `output_format` decides the file format of the output and can be `se` (default, the SE saved as RDS), `tsv` or `csv`.

## `metaproviz/internalstandard`

IS are added in the same amount to every sample. Their values should therefore not depend on the sample, and a high CV points to technical noise, for example instrument drift. Both the overall CV and the CV per condition are calculated, and the IS are removed from the output. The input data must not contain pool or QC samples.

IS are recognised by their feature names in this order: the list `--internal_standards`, the pattern `--pattern`, and the default naming conventions (`-dN`, `-13CN`, `ISTD`, `IS_` prefix and `_IS` suffix, not case sensitive). If no IS is found or every feature is an IS, the run stops. If the IS have no usable CV (all not detected, one value only or identical values), the run continues and the report and the log say which case applies.

| Parameter              | Where       | Default                                                                        | Meaning                                                                                                                                                                                              |
| ---------------------- | ----------- | ------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `output_format`        | value input | `se`                                                                           | `se`, `tsv` or `csv`.                                                                                                                                                                                |
| `--internal_standards` | `ext.args`  | none                                                                           | Comma-separated list of IS feature names. A name that contains a comma goes in double quotes.                                                                                                        |
| `--pattern`            | `ext.args`  | none                                                                           | Regular expression for the IS feature names. Used if the list is not set or finds no IS.                                                                                                             |
| `--cutoff_cv`          | `ext.args`  | `30`                                                                           | CV cutoff in percent. IS above it are flagged `HighVar`.                                                                                                                                             |
| `--conditions_col`     | `ext.args`  | `Conditions`                                                                   | Column of the sample metadata with the conditions. If it exists, a CV per condition is calculated.                                                                                                   |
| `--id_columns`         | `ext.args`  | auto-detected columns with names such as HMDB, KEGG, ChEBI, PubChem or pathway | Comma-separated columns of the feature metadata with database IDs. An IS found by the default naming conventions that has an entry there is flagged in the log and the report. The IS is still used. |

| Output              | File                                                          | Content                                                                      |
| ------------------- | ------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| `se`                | `<prefix>.is_filtered.rds`                                    | SE without the IS, written if `output_format` is `se` or empty.              |
| `filtered_matrices` | `<prefix>.is_filtered.{data,feature,sample}_matrix.{tsv,csv}` | The same data as three tables, written if `output_format` is `tsv` or `csv`. |
| `cv_table`          | `<prefix>.cv.tsv`                                             | One row per IS with `CV`, `HighVar` and `MissingValuepercentage`.            |
| `condition_cv`      | `<prefix>.condition_cv.tsv`                                   | One row per IS and one column per condition with the CV.                     |
| `plots`             | `<prefix>.plots.rds`                                          | ggplot objects, one box plot per IS and the plots of MetaProViz.             |
| `report`            | `<prefix>.report.html`                                        | HTML report.                                                                 |
| `log`               | `<prefix>.log`                                                | INFO and WARN messages.                                                      |

## Containers

All MetaProViz modules use the Bioconda package `bioconductor-metaproviz` (version 4.99.0).
