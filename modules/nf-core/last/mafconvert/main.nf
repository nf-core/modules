process LAST_MAFCONVERT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b03259f4457e393e47dfd87ea744afea462bd8614b14867e6b3640ae760f41f/data'
        : 'community.wave.seqera.io/library/last_samtools:a6d74d4fe63f646a'}"

    input:
    tuple val(meta), path(maf), val(format)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)
    tuple val(meta5), path(dict)

    output:
    tuple val(meta), path("*.axt.gz"),             optional:true, emit: axt_gz
    tuple val(meta), path("*.bam"),                optional:true, emit: bam
    tuple val(meta), path("*.bed.gz"),             optional:true, emit: bed_gz
    tuple val(meta), path("*.blast.gz"),           optional:true, emit: blast_gz
    tuple val(meta), path("*.blasttab.gz"),        optional:true, emit: blasttab_gz
    tuple val(meta), path("*.chain.gz"),           optional:true, emit: chain_gz
    tuple val(meta), path("*.cram"), path(fasta),  optional:true, emit: cram
    tuple val(meta), path("*.gff.gz"),             optional:true, emit: gff_gz
    tuple val(meta), path("*.html.gz"),            optional:true, emit: html_gz
    tuple val(meta), path("*.psl.gz"),             optional:true, emit: psl_gz
    tuple val(meta), path("*.sam.gz"),             optional:true, emit: sam_gz
    tuple val(meta), path("*.tab.gz"),             optional:true, emit: tab_gz
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    if [ -f "$dict" ]; then
        DICT_ARGS="-f ${dict}"
    else
        DICT_ARGS="-d"
    fi

    case $format in
        sam)
            maf-convert $args \$DICT_ARGS sam $maf -r 'ID:${meta.id} SM:${meta.id}' |
                samtools sort -O sam |
                gzip --no-name > ${prefix}.sam.gz
            ;;
        bam)
            maf-convert $args \$DICT_ARGS sam $maf -r 'ID:${meta.id} SM:${meta.id}' |
                samtools sort -O bam  -o ${prefix}.bam
            ;;
        cram)
            # Note 1: CRAM output is not supported if the genome is compressed with something else than bgzip.
            # Note 2: --reference is not needed because the path to the genome files is in the UR field in ${fasta}.dict
            # Note 3: To prevent relative reference path be replaced with absolute path, we disable cache and EBI querying.
            # This will not be needed in after htslib > 1.21 is released, see https://github.com/samtools/htslib/pull/1881
            export REF_CACHE='.'
            export REF_PATH='.'
            maf-convert $args \$DICT_ARGS sam $maf -r 'ID:${meta.id} SM:${meta.id}' |
                samtools sort -O cram -o ${prefix}.cram
            ;;
        *)
            maf-convert $args $format $maf |
                gzip --no-name > ${prefix}.${format}.gz
            ;;
    esac

    # maf-convert has no --version option but lastdb (part of the same package) has.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version 2>&1 | sed 's/lastdb //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    case $format in
        bam)
            touch ${prefix}.${format}
            ;;
        cram)
            touch ${prefix}.${format}
            ;;
        *)
            echo stub | gzip --no-name > ${prefix}.${format}.gz
            ;;
    esac

    # maf-convert has no --version option but lastdb (part of the same package) has.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version 2>&1 | sed 's/lastdb //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
