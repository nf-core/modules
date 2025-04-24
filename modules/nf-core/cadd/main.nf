process CADD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/biocontainers/cadd-scripts-with-envs:1.6.post1_cv1'

    containerOptions {
        if (prescored_dir) {
            ['singularity', 'apptainer'].contains(workflow.containerEngine) ?
                "-B ${annotation_dir}:/opt/CADD-scripts-1.6.post1/data/annotations -B ${prescored_dir}:/opt/CADD-scripts-1.6.post1/data/prescored" :
                "-v ${annotation_dir}:/opt/CADD-scripts-1.6.post1/data/annotations -v ${prescored_dir}:/opt/CADD-scripts-1.6.post1/data/prescored"
        } else {
            ['singularity', 'apptainer'].contains(workflow.containerEngine) ?
                "-B ${annotation_dir}:/opt/CADD-scripts-1.6.post1/data/annotations" :
                "-v ${annotation_dir}:/opt/CADD-scripts-1.6.post1/data/annotations"
        }
    }

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(annotation_dir)
    tuple val(meta3), path(prescored_dir)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.6.post1"
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    export XDG_CACHE_HOME=\$PWD/snakemake_cache
    mkdir -p \$XDG_CACHE_HOME

    cadd.sh \\
        -o ${prefix}.tsv.gz \\
        ${args} \\
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cadd: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.6.post1"
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    echo "" | gzip > ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cadd: ${VERSION}
    END_VERSIONS
    """
}
