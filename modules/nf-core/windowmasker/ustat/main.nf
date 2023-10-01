process WINDOWMASKER_USTAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::blast=2.14.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.0--h7d5a4b4_1':
        'biocontainers/blast:2.14.0--h7d5a4b4_1' }"

    input:
    tuple val(meta) , path(counts)
    tuple val(meta2), path(ref)

    output:
    tuple val(meta), path("${output}")  , emit: intervals
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    =   task.ext.args         ?: ""
    def prefix  =   task.ext.prefix       ?: "${meta.id}"
    def outfmt  =   args.contains('-outfmt fasta')                ? 'fasta'               :
                    args.contains('-outfmt maskinfo_asn1_bin')    ? 'maskinfo_asn1_bin'   :
                    args.contains('-outfmt maskinfo_asn1_text')   ? 'maskinfo_asn1_text'  :
                    args.contains('-outfmt maskinfo_xml')         ? 'maskinfo_xml'        :
                    args.contains('-outfmt seqloc_asn1_bin')      ? 'seqloc_asn1_bin'     :
                    args.contains('-outfmt seqloc_asn1_text')     ? 'seqloc_asn1_text'    :
                    args.contains('-outfmt seqloc_xml')           ? 'seqloc_xml'          :
                    'interval'

    output  = "${prefix}.${outfmt}"

    """
    windowmasker -ustat \\
        ${counts} \\
        $args \\
        -in ${ref} \\
        -out ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args    =   task.ext.args         ?: ""
    def prefix  =   task.ext.prefix       ?: "${meta.id}"
    def outfmt  =   args.contains('-outfmt fasta')                ? 'fasta'               :
                    args.contains('-outfmt maskinfo_asn1_bin')    ? 'maskinfo_asn1_bin'   :
                    args.contains('-outfmt maskinfo_asn1_text')   ? 'maskinfo_asn1_text'  :
                    args.contains('-outfmt maskinfo_xml')         ? 'maskinfo_xml'        :
                    args.contains('-outfmt seqloc_asn1_bin')      ? 'seqloc_asn1_bin'     :
                    args.contains('-outfmt seqloc_asn1_text')     ? 'seqloc_asn1_text'    :
                    args.contains('-outfmt seqloc_xml')           ? 'seqloc_xml'          :
                    'interval'

    output  = "${prefix}.${outfmt}"
    """
    touch ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
    END_VERSIONS
    """
}
