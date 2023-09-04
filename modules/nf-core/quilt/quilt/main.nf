process QUILT_QUILT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::r-quilt=1.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-quilt:1.0.4--r43h06b5641_2':
        'biocontainers/r-quilt:1.0.4--r43h06b5641_2' }"

    input:
    tuple val(meta), path(bams), path(bais), path(bamlist), val(chr), val(regions_start), val(regions_end), val(buffer), val(ngen), path(reference_haplotype_file), path(reference_legend_file), path(genetic_map_file)
    tuple val(meta2), path(posfile), path(phasefile)
    val seed
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"),              emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"),          emit: tbi,   optional:true
    tuple val(meta), path("RData", type: "dir"),    emit: rdata, optional:true
    tuple val(meta), path("plots", type: "dir"),    emit: plots, optional:true
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def extensions  = bams.collect { it.extension }
    def extension   = extensions.flatten().unique()
    def list_param  =   extension == ["bam"]  ? "bamlist=\"${bamlist}\""   :
                        extension == ["cram"] ? "cramlist=\"${bamlist}\", reference=\"${fasta}\"" :
                        ""

    """


    R -e 'library("QUILT") ; QUILT(
        chr="$chr",
        $list_param,
        posfile="$posfile",
        phasefile="$phasefile",
        reference_haplotype_file="$reference_haplotype_file",
        genetic_map_file="$genetic_map_file",
        nGen=$ngen,
        regionStart=$regions_start,
        regionEnd=$regions_end,
        buffer=$buffer,
        outputdir="." $args,
        reference_legend_file="$reference_legend_file",
        nCores=$task.cpus,  # Using more than 1 core may lead to different md5sum
        seed=$seed
        )'


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-quilt: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"QUILT\\")))")
    END_VERSIONS
    """
}
