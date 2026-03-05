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
    tuple val(meta), path("*.{axt.gz,bam,bed.gz,blast.gz,blasttab.gz,chain.gz,cram,gff.gz,html.gz,psl.gz,sam.gz,tab.gz}"), emit: alignment
    // last-dotplot has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('last'), eval("lastal --version | sed 's/lastal //'"), emit: versions_last, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    dict2gff3() { awk '
        BEGIN {
            print "##gff-version 3"
        }
        \$1 == "@SQ" {
            seq_name   = ""
            seq_length = ""
            for (i = 1; i <= NF; i++) {
                if      (\$i ~ /^SN:/) seq_name   = substr(\$i, 4)
                else if (\$i ~ /^LN:/) seq_length = substr(\$i, 4)
            }
            if (seq_name != "" && seq_length != "") {
                printf "##sequence-region %s 1 %s\\n", seq_name, seq_length
            }
        }' "\${1:-/dev/stdin}"
    }

    if [ -f "$dict" ]; then
        DICT_ARGS="-f ${dict}"
        [ "$format" = "gff" ] && dict2gff3 ${dict}          > "${prefix}.head.gff"
    else
        DICT_ARGS="-d"
        [ "$format" = "gff" ] && printf "##gff-version 3\\n" > "${prefix}.head.gff"
    fi

    case $format in
        gff)
            cat "${prefix}.head.gff" <(maf-convert $args -n gff $maf) |
                gzip --no-name > ${prefix}.gff.gz
            ;;
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
            # Note 4: CRAM version 3.0 is enforced until htsjdk, and therefore nf-test, supports 3.1
            maf-convert $args \$DICT_ARGS sam $maf -r 'ID:${meta.id} SM:${meta.id}' |
                samtools sort -O cram,version=3.0 -o ${prefix}.cram
            ;;
        *)
            maf-convert $args $format $maf |
                gzip --no-name > ${prefix}.${format}.gz
            ;;
    esac
    """

    stub:
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
    """
}
