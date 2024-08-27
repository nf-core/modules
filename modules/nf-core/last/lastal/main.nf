process LAST_LASTAL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1571--h43eeafb_0' :
        'biocontainers/last:1571--h43eeafb_0' }"

    input:
    tuple val(meta), path(fastx), path (param_file)
    path index

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    tuple val(meta), path("*.tsv")   , emit: multiqc
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trained_params = param_file ? "-p ${param_file}"  : ''
    """
    INDEX_NAME=\$(basename \$(ls $index/*.des) .des)
    set -o pipefail

    function calculate_psl_metrics() {
        awk 'BEGIN {
            FS="\t";  # Set field separator as tab
            totalMatches = 0;
            totalAlignmentLength = 0;
            print "Sample\tTotalAlignmentLength\tPercentSimilarity";  # Header for MultiQC
        }
        {
            totalMatches += \$1 + \$3;  # Sum matches and repMatches
            totalAlignmentLength += \$1 + \$2 + \$3 + \$6 + \$8;  # Sum matches, misMatches, repMatches, qBaseInsert, and tBaseInsert
        }
        END {
            percentSimilarity = (totalAlignmentLength > 0) ? (totalMatches / totalAlignmentLength * 100) : 0;
            print "$meta.id" "\t" totalAlignmentLength "\t" percentSimilarity;  # Data in TSV format
        }'
    }

    lastal \\
        -P $task.cpus \\
        $trained_params \\
        $args \\
        ${index}/\$INDEX_NAME \\
        $fastx |
        tee >(gzip --no-name  > ${prefix}.maf.gz) |
        maf-convert psl |
        calculate_psl_metrics > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastal --version 2>&1 | sed 's/lastal //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trained_params = param_file ? "-p ${param_file}"  : ''
    """
    INDEX_NAME=STUB
    echo stub | gzip --no-name > ${prefix}.\$INDEX_NAME.maf.gz
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastal --version 2>&1 | sed 's/lastal //')
    END_VERSIONS
    """
}
