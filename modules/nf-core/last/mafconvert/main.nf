process LAST_MAFCONVERT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/8975c4b3f5acfdf065246236c3a1b8cde7a0df4cacff72f9ca2670fc7626726e/data'
        : 'community.wave.seqera.io/library/bcftools_last_samtools:f6cbdfb2676292ee'}"

    input:
    tuple val(meta), path(maf), val(format)
    tuple val(meta2), path(fasta), path(fai), path(gzi), path(sizes), path(dict) // see subworkflows/nf-core/fasta_bgzip_index_dict_samtools

    output:
    tuple val(meta), path("*.{axt.gz,bam,bcf,bed.gz,blast.gz,blasttab.gz,chain.gz,cram,gff.gz,html.gz,psl.gz,sam.gz,tab.gz}"), emit: alignment
    tuple val(meta), path("*.{bai,crai,csi}"), emit: index, optional: true
    tuple val(meta), path("*.stats"),          emit: stats, optional: true
    // last-dotplot has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('last'), eval("lastal --version | sed 's/lastal //'"), emit: versions_last, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''   // maf-convert
    def args2  = task.ext.args2  ?: ''   // samtools sort
    def args3  = task.ext.args3  ?: ''   // bcftools mpileup
    def args4  = task.ext.args4  ?: ''   // bcftools call
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    if [ -f "$dict" ]; then
        DICT_ARGS="-f ${dict}"
        if [ "$format" = "cram" ]; then
            REF_CRAM=\$(grep '^@SQ' $dict | sed -n 's/.*UR:\\([^ \\t]*\\).*/\\1/p' | uniq)
            if [ -r \$REF_CRAM ]; then
                REF_ARGS=''
            else
                REF_ARGS="--reference $fasta"
            fi
        fi
    else
        DICT_ARGS="-d"
    fi

    case $format in
        gff)
            {
                echo "##gff-version 3"
                [ -f "$sizes" ] && awk '{ printf "##sequence-region %s 1 %s\\n", \$1, \$2 }' $sizes
                maf-convert $args -n gff $maf
            } | gzip --no-name > ${prefix}.gff.gz
            ;;
        sam)
            maf-convert $args \$DICT_ARGS sam $maf -r 'ID:${meta.id} SM:${meta.id}' |
                samtools sort -O sam |
                gzip --no-name > ${prefix}.sam.gz
            ;;
        bam)
            maf-convert $args \$DICT_ARGS sam $maf -r 'ID:${meta.id} SM:${meta.id}' |
                samtools sort $args2 -O bam  -o ${prefix}.bam
            ;;
        cram)
            # Note 1: CRAM output is not supported if the genome is compressed with something else than bgzip.
            maf-convert $args \$DICT_ARGS sam $maf -r 'ID:${meta.id} SM:${meta.id}' |
                samtools sort $args2 -O cram \$REF_ARGS -o ${prefix}.cram
            ;;
        bcf)
            maf-convert $args \$DICT_ARGS sam $maf -r 'ID:${meta.id} SM:${meta.id}' |
                samtools sort $args2 -u | bcftools mpileup $args3 --fasta-ref $fasta -Ou - | bcftools call $args4 -Ob > ${prefix}.bcf
            bcftools stats ${prefix}.bcf > ${prefix}.stats
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
            echo "" | gzip > ${prefix}.${format}.gz
            ;;
    esac
    """
}
