process PLINK_RECODE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.ped")                    , optional:true, emit: ped
    tuple val(meta), path("*.map")                    , optional:true, emit: map
    tuple val(meta), path("*.txt")                    , optional:true, emit: txt
    tuple val(meta), path("*.raw")                    , optional:true, emit: raw
    tuple val(meta), path("*.traw")                   , optional:true, emit: traw
    tuple val(meta), path("*.beagle.dat")             , optional:true, emit: beagledat
    tuple val(meta), path("*.chr-*.dat")              , optional:true, emit: chrdat
    tuple val(meta), path(".*chr-*.map")              , optional:true, emit: chrmap
    tuple val(meta), path("*.recode.geno.txt")        , optional:true, emit: geno
    tuple val(meta), path("*.recode.pheno.txt")       , optional:true, emit: pheno
    tuple val(meta), path("*.recode.pos.txt")         , optional:true, emit: pos
    tuple val(meta), path("*.recode.phase.inp")       , optional:true, emit: phase
    tuple val(meta), path("*.info")                   , optional:true, emit: info
    tuple val(meta), path("*.lgen")                   , optional:true, emit: lgen
    tuple val(meta), path("*.list")                   , optional:true, emit: list
    tuple val(meta), path("*.gen")                    , optional:true, emit: gen
    tuple val(meta), path("*.gen.gz")                 , optional:true, emit: gengz
    tuple val(meta), path("*.sample")                 , optional:true, emit: sample
    tuple val(meta), path("*.rlist")                  , optional:true, emit: rlist
    tuple val(meta), path("*.strct_in")               , optional:true, emit: strctin
    tuple val(meta), path("*.tped")                   , optional:true, emit: tped
    tuple val(meta), path("*.tfam")                   , optional:true, emit: tfam
    tuple val(meta), path("*.vcf")                    , optional:true, emit: vcf
    tuple val(meta), path("*.vcf.gz")                 , optional:true, emit: vcfgz
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plink \\
        --bed ${bed}  \\
        --bim ${bim}  \\
        --fam ${fam}  \\
        --threads ${task.cpus} \\
        --recode \\
        ${args} \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.ped
    touch ${prefix}.map

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """
}
