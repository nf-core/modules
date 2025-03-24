process MINIMAC4_IMPUTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimac4:4.1.6--hcb620b3_1':
        'biocontainers/minimac4:4.1.6--hcb620b3_1' }"

    input:
    tuple val(meta), path(target_vcf), path(target_index), path(ref_msav), path(sites_vcf), path(sites_index), path(map)

    output:
    tuple val(meta), path("*.{bcf,sav,vcf.gz,vcf,ubcf,usav}"), emit: vcf
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-format bcf")    || args.contains("-O bcf")    ? "bcf"    :
                    args.contains("--output-format sav")    || args.contains("-O sav")    ? "sav"    :
                    args.contains("--output-format vcf.gz") || args.contains("-O vcf.gz") ? "vcf.gz" :
                    args.contains("--output-format vcf")    || args.contains("-O vcf")    ? "vcf"    :
                    args.contains("--output-format ubcf")   || args.contains("-O ubcf")   ? "ubcf"   :
                    args.contains("--output-format usav")   || args.contains("-O usav")   ? "usav"   :
                    "vcf.gz"
    def sites_cmd = sites_vcf ? "--sites $sites_vcf" : ""
    def map_cmd   = map       ? "--map $map"         : ""
    """
    minimac4 \\
        $ref_msav \\
        $target_vcf \\
        $args \\
        $sites_cmd \\
        $map_cmd \\
        --threads $task.cpus \\
        -o ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: \$(minimac4 --version |& sed '1!d ; s/minimac v//')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-format bcf")    || args.contains("-O bcf")    ? "bcf"    :
                    args.contains("--output-format sav")    || args.contains("-O sav")    ? "sav"    :
                    args.contains("--output-format vcf.gz") || args.contains("-O vcf.gz") ? "vcf.gz" :
                    args.contains("--output-format vcf")    || args.contains("-O vcf")    ? "vcf"    :
                    args.contains("--output-format ubcf")   || args.contains("-O ubcf")   ? "ubcf"   :
                    args.contains("--output-format usav")   || args.contains("-O usav")   ? "usav"   :
                    "vcf.gz"
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    """
    ${create_cmd} ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimac4: \$(minimac4 --version |& sed '1!d ; s/minimac v//')
    END_VERSIONS
    """
}
