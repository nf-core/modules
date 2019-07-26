/*
* Description:
*     Run FastQC on sequenced reads
* Keywords:
*     read qc
*     adapter
* Tools:
*     FastQC:
*         homepage: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
*         documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/
*         description: FastQC gives general quality metrics about your reads.
*                      It provides information about the quality score distribution
*                      across your reads, the per base sequence content (%A/C/G/T).
*                      You get information about adapter contamination and other
*                      overrepresented sequences.
*/
process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(sample_id), file(reads)

    output:
    file "*_fastqc.{zip,html}"

    script:
    """
    fastqc -q $reads

    fastqc --version &> fastqc.version.txt
    """
}
