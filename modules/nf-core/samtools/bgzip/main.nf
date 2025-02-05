process SAMTOOLS_BGZIP {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/61/61be440cd54169fcefe5303238847dd466728453a130f4fbc5abd68b514f0b09/data' :
        'community.wave.seqera.io/library/samtools_file:f35ca613dd912bed' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.bgzip.fa.gz") , emit: fa
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    COMPRESS_TYPE=\$(file -L -b $fasta)

    case "\$COMPRESS_TYPE" in
        "gzip compressed data, extra field"*)
            # A well-behaved find installation should report:
            # Blocked GNU Zip Format (BGZF; gzip compatible)
            # But the one in anaconda does nt…
            # Assuming the "extra field" implies BGZF, do nothing.
            ln -s $fasta ${prefix}.bgzip.fa.gz
            ;;
        gzip*)
            # Recompress non-BGZF gzipped files
            zcat $fasta |
                bgzip \\
                    $args \\
                    --threads ${task.cpu} \\
                    > ${prefix}.bgzip.fa.gz
            ;;
        *)
            # Compress
            bgzip \\
                $args \\
                --threads ${task.cpu} \\
                --stdout \\
                $fasta \\
                > ${prefix}.bgzip.fa.gz
            ;;
    esac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        file: \$(echo \$(file --version 2>&1 | head -n1 | sed 's/file-//'))
        samtools: \$(echo \$(samtools --version 2>&1 | head -n1 | sed 's/^.*samtools //; s/Using.*\$//'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | bgzip > ${prefix}.bgzip.fa.gz

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        file: \$(echo \$(file --version 2>&1 | head -n1 | sed 's/file-//'))
        samtools: \$(echo \$(samtools --version 2>&1 | head -n1 | sed 's/^.*samtools //; s/Using.*\$//'))
    END_VERSIONS
    """
}
