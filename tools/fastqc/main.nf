process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    container: 'quay.io/biocontainers/fastqc:0.11.8--2'

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
