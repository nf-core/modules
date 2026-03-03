process QUILT_QUILT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/r-quilt:1.0.5--r43h06b5641_0'
        : 'biocontainers/r-quilt:1.0.5--r43h06b5641_0'}"

    input:
    tuple val(meta), path(bams), path(bais), path(bamlist), path(samplename), path(reference_haplotype_file), path(reference_legend_file), path(posfile), path(phasefile), path(genfile), val(chr), val(regions_start), val(regions_end), val(ngen), val(buffer), path(genetic_map)
    tuple val(meta2), path(fasta), path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz")          , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")      , emit: tbi  , optional: true
    tuple val(meta), path("RData", type: "dir"), emit: rdata, optional: true
    tuple val(meta), path("plots", type: "dir"), emit: plots, optional: true
    tuple val("${task.process}"), val('r-quilt'), eval('Rscript -e "cat(as.character(packageVersion(\'QUILT\')))"'), topic: versions, emit: versions_r_quilt
    tuple val("${task.process}"), val('r-base'), eval('R --version | sed "1!d; s/.*version //; s/ .*//"'), topic: versions, emit: versions_r_base

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"

    def extensions   = bams.collect { path -> path.extension }
    def extension    = extensions.flatten().unique()
    def list_command = extension == ["bam"]
        ? "--bamlist="
        : extension == ["cram"] ? "--reference=${fasta} --cramlist=" : ""

    def genetic_map_command = genetic_map   ? "--genetic_map_file=${genetic_map}" : ""
    def posfile_command     = posfile       ? "--posfile=${posfile}"              : ""
    def phasefile_command   = phasefile     ? "--phasefile=${phasefile}"          : ""
    def samplename_command  = samplename    ? "--sampleNames_file=${samplename}"  : ""
    def start_command       = regions_start ? "--regionStart=${regions_start}"    : ""
    def end_command         = regions_end   ? "--regionEnd=${regions_end}"        : ""
    def buffer_command      = buffer        ? "--buffer=${buffer}"                : ""

    if (!(args ==~ /.*--seed.*/)) {
        args += " --seed=1"
    }

    """
    if [ -n "${bamlist}" ] ;
    then
        BAM_LIST="${bamlist}"
    else
        printf "%s\\n" ${bams} | tr -d '[],' > all_files.txt
        BAM_LIST="all_files.txt"
    fi

    QUILT.R \\
        ${list_command}\$BAM_LIST \\
        ${genetic_map_command} \\
        ${posfile_command} \\
        ${phasefile_command} \\
        ${samplename_command} \\
        --chr=${chr} \\
        ${start_command} \\
        ${end_command} \\
        ${buffer_command} \\
        --nGen=${ngen} \\
        --nCores=${task.cpus} \\
        --outputdir="." \\
        --reference_haplotype_file=${reference_haplotype_file} \\
        --reference_legend_file=${reference_legend_file} \\
        --output_filename=${prefix}.${suffix} \\
        ${args}
    """

    stub:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def suffix        = task.ext.suffix ?: "vcf.gz"
    def create_cmd    = suffix.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def make_plots    = args.contains("--make_plots=TRUE")
    def save_ref      = args.contains("--save_prepared_reference=TRUE")
    def nGibbsSamples = args.contains("--nGibbsSamples=") ? args.split("--nGibbsSamples=")[1].split(" ")[0] : 7
    def n_seek_its    = args.contains("--n_seek_its=")    ? args.split("--n_seek_its=")[1].split(" ")[0]    : 3

    """
    ${create_cmd} ${prefix}.${suffix}
    touch ${prefix}.${suffix}.tbi
    if [ "${save_ref}" == true ]
    then
        mkdir -p RData
        touch "RData/QUILT_prepared_reference.${chr}.${regions_start}.${regions_end}.RData"
    fi
    if [ "${make_plots}" == true ]
    then
        mkdir -p plots
        for nGibbs in {0..${nGibbsSamples}}
        do
            touch "plots/haps.${prefix}.${chr}.${regions_start}.${regions_end}_igs.\$((nGibbs+1)).0.truth.png"
            for its in {1..${n_seek_its}}
            do
                touch "plots/haps.${prefix}.${chr}.${regions_start}.${regions_end}_igs.\$((nGibbs+1)).it\$its.gibbs.png"
            done
        done
    fi
    """
}
