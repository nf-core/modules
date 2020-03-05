process tcoffee {
    tag {fasta}
    // publishDir "${params.outdir}/tcoffee" //not working
    publishDir "${baseDir}/results/tcoffee", mode: 'copy'
    container 'quay.io/biocontainers/t_coffee:11.0.8--py27pl5.22.0_5'

    input:
    path(fasta)

    output:
    path("${fasta}.aln")

    script:
    """
    t_coffee -seq $fasta -outfile ${fasta}.aln
    """
}

