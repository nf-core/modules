process LAST_SPLIT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db0b5de918238f07ec1ca668be942397da85e26aa582f8927ac37c70896303cf/data'
        : 'community.wave.seqera.io/library/last:1608--f41c047f7dc37e30'}"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    tuple val(meta), path("*.tsv")   , emit: multiqc
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
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

    zcat < $maf |
        last-split $args |
        tee >(gzip --no-name  > ${prefix}.maf.gz) |
        maf-convert psl |
        calculate_psl_metrics > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(last-split --version 2>&1 | sed 's/last-split //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$maf" == "${prefix}.maf.gz" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo stub | gzip --no-name > ${prefix}.maf.gz
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(last-split --version 2>&1 | sed 's/last-split //')
    END_VERSIONS
    """
}
