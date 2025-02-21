process QUILT_QUILT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-quilt:1.0.5--r43h06b5641_0':
        'biocontainers/r-quilt:1.0.5--r43h06b5641_0' }"

    input:
    tuple val(meta), path(bams), path(bais), path(bamlist), path(samplename), path(reference_haplotype_file), path(reference_legend_file), val(chr), val(regions_start), val(regions_end), val(ngen), val(buffer), path(genetic_map_file)
    tuple val(meta2), path(posfile), path(phasefile), path(genfile)
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
    def args                        =   task.ext.args ?: ''
    def _prefix                     =   task.ext.prefix ?: "${meta.id}"
    def extensions                  =   bams.collect { it.extension }
    def extension                   =   extensions.flatten().unique()
    def list_command                =   extension == ["bam"]  ? "--bamlist="                       :
                                        extension == ["cram"] ? "--reference=${fasta} --cramlist=" : ""
    def genetic_map_file_command    =   genetic_map_file      ? "--genetic_map_file=${genetic_map_file}"     : ""
    def posfile_command             =   posfile               ? "--posfile=${posfile}"                       : ""
    def phasefile_command           =   phasefile             ? "--phasefile=${phasefile}"                   : ""
    def samplename_command          =   samplename            ? "--sampleNames_file=${samplename}"           : ""
    if (!(args ==~ /.*--seed.*/)) {args += " --seed=1"}

    """
    if [ -n "$bamlist" ] ;
    then
        BAM_LIST="$bamlist"
    else
        printf "%s\\n" $bams | tr -d '[],' > all_files.txt
        BAM_LIST="all_files.txt"
    fi

    QUILT.R \\
        ${list_command}\$BAM_LIST \\
        $genetic_map_file_command \\
        $posfile_command \\
        $phasefile_command \\
        $samplename_command \\
        --chr=$chr \\
        --regionStart=$regions_start \\
        --regionEnd=$regions_end \\
        --nGen=$ngen \\
        --buffer=$buffer \\
        --nCores=$task.cpus \\
        --outputdir="." \\
        --reference_haplotype_file=$reference_haplotype_file \\
        --reference_legend_file=$reference_legend_file \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-quilt: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"QUILT\\")))")
    END_VERSIONS
    """
}
