process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda 'bioconda::clair3=1.0.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.0.0--py39hb9dc472_1' :
        'quay.io/biocontainers/clair3:1.0.0--py39hb9dc472_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}/*")        , emit: vcf
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    /usr/local/bin/run_clair3.sh \
    --bam_fn=$bam \
    --ref_fn=$fasta \
    --threads=$task.cpus \
    --output="${meta.id}" \
    $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$( /usr/local/bin/run_clair3.sh --version | sed 's/Clair3 v//' )
    END_VERSIONS
    """
}
