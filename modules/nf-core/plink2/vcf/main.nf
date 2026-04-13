process PLINK2_VCF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a42bdfda9d3e7e247cfa2e8cda32a3e40f8fa8c1a5e9bce8b8b7b4c8fd6d3f49/data' :
        'community.wave.seqera.io/library/plink2:2.0a6.9--b1e8b16e8fc23b39' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.pgen")    , emit: pgen
    tuple val(meta), path("*.psam")    , emit: psam
    tuple val(meta), path("*.pvar")    , emit: pvar
    tuple val(meta), path("*.pvar.zst"), emit: pvar_zst, optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --make-pgen \\
        --set-all-var-ids '@:#:\$r:\$a' \\
        --new-id-max-allele-len 10 missing \\
        --rm-dup force-first \\
        $args \\
        --vcf $vcf \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
