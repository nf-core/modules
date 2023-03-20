process PEPPER_MARGIN_DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    if (params.deepvariant_gpu) {
        container 'docker.io/kishwars/pepper_deepvariant:r0.8-gpu'
    } else {
        container 'docker.io/kishwars/pepper_deepvariant:r0.8'
    }

    input:
    tuple val(meta), path(input), path(index), val(intervals)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*vcf.gz")    ,  emit: vcf
    tuple val(meta), path("*vcf.gz.tbi"),  emit: tbi
    path "versions.yml"                 ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ""
    def gpu     = params.deepvariant_gpu ? "-g" : ""
    prefix      = task.ext.prefix ?: "${meta.id}"
    //def regions = intervals ? "--regions ${intervals}" : ""
    //def gvcf    = params.make_gvcf ? "--gvcf" : ""

    """
    mkdir -p "${prefix}"
    run_pepper_margin_deepvariant call_variant \\
        -b "${input}" \\
        -f "${fasta}" \\
        -o "." \\
        -p "${prefix}" \\
        -t ${task.cpus} \\
        $gpu \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepper_margin_deepvariant: \$(run_pepper_margin_deepvariant --version | sed 's/VERSION: //g')
    END_VERSIONS
    """
}
