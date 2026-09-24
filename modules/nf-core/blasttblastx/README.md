# blasttblastx

This module runs **NCBI BLAST+ `tblastx`**, comparing six-frame translations of nucleotide sequences in a query FASTA file against a six-frame translated nucleotide database.

> `tblastx` is used to detect coding regions and conserved domains between distantly related nucleotide sequences where protein-level homology is more detectable than nucleotide-level.

---

## Usage

This module expects the following inputs:

- A nucleotide **query FASTA** file (`params.query`)
- A **BLAST nucleotide database** directory (`params.db`) built using `makeblastdb -dbtype nucl`

The output is:

- A BLAST **tabular results file** (`tblastx_results.txt`) in outfmt 6

---

## Example

```bash
nextflow run main.nf -c tests/config/nextflow.config