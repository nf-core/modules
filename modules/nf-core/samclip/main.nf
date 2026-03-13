process SAMCLIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/samclip_samtools:7af2916e4ae6f461'
        : 'community.wave.seqera.io/library/samclip_samtools:00cc7aefd75be672'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(reference), path(reference_index)

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: reads
    tuple val("${task.process}"), val('samclip'), eval("samclip --version | sed 's/^.*samclip //g'"), topic: versions, emit: versions_samclip
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''   // samclip args
    def args2      = task.ext.args2  ?: ''   // samtools sort (first, name sort) args
    def args3      = task.ext.args3  ?: ''   // samtools fixmate args
    def args4      = task.ext.args4  ?: ''   // samtools sort (second, coordinate sort) args
    def prefix     = task.ext.prefix ?: "${meta.id}_samclip"
    def extension  = args4.contains("--output-fmt cram") ? "cram" :
                     args4.contains("-O cram")           ? "cram" :
                     args4.contains("-O CRAM")           ? "cram" :
                     "bam"
    def reference_arg = extension == "cram" ? "--reference ${reference}" : ""
    def is_compressed = reference.getName().endsWith(".gz")
    def ref_filename  = reference.getName().replaceAll(/\.gz$/, "")

    if ("${bam}" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    # decompress reference if gzipped
    if [ "${is_compressed}" = "true" ]; then
        gzip -c -d ${reference} > ${ref_filename}
    fi

    samtools view -h --output-fmt sam ${bam} | \\
    samclip ${args} --ref ${ref_filename} | \\
    samtools sort -n -O SAM ${args2} | \\
    samtools fixmate -m ${args3} - - | \\
    samtools sort ${args4} ${reference_arg} -O ${extension.toUpperCase()} -o ${prefix}.${extension}

    # clean up decompressed reference
    if [ "${is_compressed}" = "true" ]; then
        rm -f ${ref_filename}
    fi

    """

    stub:
    def args4     = task.ext.args4  ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}_samclip"
    def extension = args4.contains("--output-fmt cram") ? "cram" :
                    args4.contains("-O cram")           ? "cram" :
                    args4.contains("-O CRAM")           ? "cram" :
                    "bam"

    if ("${bam}" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.${extension}

    """
}
