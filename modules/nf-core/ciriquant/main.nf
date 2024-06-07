process CIRIQUANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ciriquant=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciriquant:1.1.2--pyhdfd78af_2' :
        'biocontainers/ciriquant:1.1.2--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(gtf)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(bwa)
    tuple val(meta5), path(hisat2)

    output:
    tuple val(meta), path("${prefix}/${prefix}.gtf"), emit: gtf
    path("${prefix}")                               , emit: results
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.0'
    """
    BWA=`which bwa`
    HISAT2=`which hisat2`
    STRINGTIE=`which stringtie`
    SAMTOOLS=`which samtools`

    BWA_FILE=`ls ${bwa}/*.bwt`
    BWA_PREFIX=`basename \$BWA_FILE .bwt`

    HISAT2_FILE=`ls ${hisat2}/*.1.ht2`
    HISAT2_PREFIX=`basename \$HISAT2_FILE .1.ht2`

    printf "name: ciriquant\\ntools:\\n  bwa: \$BWA\\n  hisat2: \$HISAT2\\n  stringtie: \$STRINGTIE\\n  samtools: \$SAMTOOLS\\n\\nreference:\\n  fasta: ${fasta}\\n  gtf: ${gtf}\\n  bwa_index: ${bwa}/\$BWA_PREFIX\\n  hisat_index: ${hisat2}/\$HISAT2_PREFIX" > config.yml

    CIRIquant \\
        -t ${task.cpus} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --config config.yml \\
        --no-gene \\
        -o ${prefix} \\
        -p ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        ciriquant: \$(echo \$(CIRIquant --version 2>&1) | sed 's/CIRIquant //g' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        stringtie: \$(stringtie --version 2>&1)
        hisat2: $VERSION
    END_VERSIONS
    """
}
