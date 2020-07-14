process shovill {

    tag { shovill }

    publishDir "${params.outdir}", pattern: '*.fasta', mode: 'copy'

    container "quay.io/biocontainers/shovill:1.0.9--0"

    input:
    tuple(sample_id, path(forward), path(reverse))

    output:
    path("${sample_id}.fasta")

    script:
    """
    shovill --R1 ${forward} --R2 ${reverse} --outdir shovill_out
    mv shovill_out/contigs.fa ${sample_id}.fasta
    """
}
