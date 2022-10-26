process HIFIASM {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::hifiasm=0.15.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hifiasm:0.15.4--h2e03b76_0' :
        'quay.io/biocontainers/hifiasm:0.15.4--h2e03b76_0' }"

    input:
    tuple val(meta), path(reads)
    path  paternal_kmer_dump
    path  maternal_kmer_dump
    val   use_parental_kmers

    output:
    tuple val(meta), path("*.r_utg.gfa")       , emit: raw_unitigs
    tuple val(meta), path("*.ec.bin")          , emit: corrected_reads
    tuple val(meta), path("*.ovlp.source.bin") , emit: source_overlaps
    tuple val(meta), path("*.ovlp.reverse.bin"), emit: reverse_overlaps
    tuple val(meta), path("*.p_utg.gfa")       , emit: processed_unitigs, optional: true
    tuple val(meta), path("*.asm.p_ctg.gfa")   , emit: primary_contigs  , optional: true
    tuple val(meta), path("*.asm.a_ctg.gfa")   , emit: alternate_contigs, optional: true
    tuple val(meta), path("*.hap1.p_ctg.gfa")  , emit: paternal_contigs , optional: true
    tuple val(meta), path("*.hap2.p_ctg.gfa")  , emit: maternal_contigs , optional: true
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (use_parental_kmers) {
        """
        hifiasm \\
            $args \\
            -o ${prefix}.asm \\
            -t $task.cpus \\
            -1 $paternal_kmer_dump \\
            -2 $maternal_kmer_dump \\
            $reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """
    } else { // Phasing with Hi-C data is not supported yet
        """
        hifiasm \\
            $args \\
            -o ${prefix}.asm \\
            -t $task.cpus \\
            $reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """
    }
}
