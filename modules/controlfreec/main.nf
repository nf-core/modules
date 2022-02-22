process CONTROLFREEC_SOMATIC {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::control-freec=11.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/control-freec:11.6--h1b792b2_1':
        'quay.io/biocontainers/control-freec:11.6--h1b792b2_1' }"

    input:
    tuple val(meta), path(mpileup_normal), path(mpileup_tumor), path(minipileup_normal), path(minipileup_tumor)
    tuple val(meta), path(cpn_normal), path(cpn_tumor), path(minipileup_normal), path(minipileup_tumor)
    path fasta
    path fasta_fai
    path known_snps
    path known_snps_tbi
    path chr_length
    path chr_directory
    path mappability
    path target_bed

    output:
    tuple val(meta), path("*.txt"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    //"General" configurations
    def bedgraphoutput              = task.ext.args?["general"]?["bedgraphoutput"]              ? "BedGraphOutput = ${task.ext.args["general"]["bedgraphoutput"]}"                              : ""
    //bedtools: not needed since pileup files are defined as input
    def breakpointthreshold         = task.ext.args?["general"]?["breakpointthreshold"]         ? "breakPointThreshold = ${task.ext.args["general"]["breakpointthreshold"]}"                    : ""
    def breakpointtype              = task.ext.args?["general"]?["breakpointtype"]              ? "breakPointType = ${task.ext.args["general"]["breakpointtype"]}"                              : ""
    def coefficientofvariation      = task.ext.args?["general"]?["coefficient"]                 ? "coefficientOfVariation = ${task.ext.args["general"]["coefficientofvariation"]}"              : ""
    def contamination               = task.ext.args?["general"]?["contamination"]               ? "contamination = ${task.ext.args["general"]["contamination"]}"                                : ""
    def contaminationadjustment     = task.ext.args?["general"]?["contaminationadjustment"]     ? "contaminationAdjustment = ${task.ext.args["general"]["contaminationadjustment"]}"            : ""
    def degree                      = task.ext.args?["general"]?["degree"]                      ? "degree = ${task.ext.args["general"]["degree"]}"                                              : ""
    def forcegccontentnormalization = task.ext.args?["general"]?["forcegccontentnormalization"] ? "forceGCcontentNormalization = ${task.ext.args["general"]["forcegccontentnormalization"]}"    : ""
    def gccontentprofile            = task.ext.args?["general"]?["gccontentprofile"]            ? "GCcontentProfile = ${task.ext.args["general"]["gccontentprofile"]}"                          : ""
    def mappability                 = mappability                                               ? "gemMappabilityFile = \${PWD}/${mappability}"                                                 : ""
    def intercept                   = task.ext.args?["general"]?["intercept"]                   ? "intercep = ${task.ext.args["general"]["intercept"]}"                                         : ""
    def mincnalength                = task.ext.args?["general"]?["mincnalength"]                ? "minCNAlength = ${task.ext.args["general"]["mincnalength"]}"                                  : ""
    def minmappabilityperwindow     = task.ext.args?["general"]?["minmappabilityperwindow"]     ? "minMappabilityPerWindow = ${task.ext.args["general"]["minmappabilityperwindow"]}"            : ""
    def minexpectedgc               = task.ext.args?["general"]?["minexpectedgc"]               ? "minExpectedGC = ${task.ext.args["general"]["minexpectedgc"]}"                                : ""
    def maxexpectedgc               = task.ext.args?["general"]?["maxexpectedgc"]               ? "maxExpectedGC = ${task.ext.args["general"]["maxexpectedgc"]}"                                : ""
    def minimalsubclonepresence     = task.ext.args?["general"]?["minimalsubclonepresence"]     ? "minimalSubclonePresence = ${task.ext.args["general"]["minimalsubclonepresence"]}"            : ""
    def noisydata                   = task.ext.args?["general"]?["noisydata"]                   ? "noisyData = ${task.ext.args["general"]["noisydata"]}"                                        : ""
    def ploidy                      = task.ext.args?["general"]?["ploidy"]                      ? "ploidy = ${task.ext.args["general"]["ploidy"]}"                                              : ""
    def printNA                     = task.ext.args?["general"]?["printNA"]                     ? "printNA = ${task.ext.args["general"]["printNA"]}"                                            : ""
    def readcountthreshold          = task.ext.args?["general"]?["readcountthreshold"]          ? "readCountThreshold = ${task.ext.args["general"]["readcountthreshold"]}"                      : ""
    def sex                         = task.ext.args?["general"]?["sex"]                         ? "sex = ${task.ext.args["general"]["sex"]}"                                                    : ""
    def step                        = task.ext.args?["general"]?["step"]                        ? "step = ${task.ext.args["general"]["step"]}"                                                  : ""
    def telocentromeric             = task.ext.args?["general"]?["telocentromeric"]             ? "telocentromeric = ${task.ext.args["general"]["telocentromeric"]} "                           : ""
    def uniquematch                 = task.ext.args?["general"]?["uniquematch"]                 ? "uniqueMatch = ${task.ext.args["general"]["uniquematch"]}"                                    : ""
    def window                      = task.ext.args?["general"]?["window"]                      ? "window = ${task.ext.args["general"]["window"]}"                                              : ""

    //"Control" configurations
    def matefile_normal             = mpileup_normal                                            ? "mateFile = \${PWD}${mpileup_normal}"                                                         : ""
    def matecopynumberfile_normal   = cpn_normal                                                ? "mateCopyNumberFile = \${PWD}${cpn_normal}"                                                   : ""
    def minipileup_normal           = minipileup_normal                                         ? "miniPileup = \${PWD}${minipileup_normal}"                                                    : ""
    def inputformat_normal          = task.ext.args?["general"]?["inputformat_normal"]          ? "inputFormat = ${task.ext.args["general"]["inputformat_normal"]}"                             : ""
    def mateorientation_normal      = task.ext.args?["general"]?["mateorientation_normal"]      ? "mateOrientation = ${task.ext.args["general"]["mateorientation_normal"]}"                     : ""

    //"Sample" configuration
    def matefile_tumor             = mpileup_tumor                                              ? "mateFile = \${PWD}${mpileup_tumor}"                                                          : ""
    def matecopynumberfile_tumor   = cpn_tumor                                                  ? "mateCopyNumberFile = \${PWD}${cpn_tumor}"                                                    : ""
    def minipileup_tumor           = minipileup_tumor                                           ? "miniPileup = \${PWD}${minipileup_tumor}"                                                     : ""
    def inputformat_tumor          = task.ext.args?["general"]?["inputformat_tumor"]            ? "inputFormat = ${task.ext.args["general"]["inputformat_tumor"]}"                              : ""
    def mateorientation_tumor      = task.ext.args?["general"]?["mateorientation_tumor"]        ? "mateOrientation = ${task.ext.args["general"]["mateorientation_tumor"]}"                      : ""


    //"BAF" configuration
    def makepileup                 =
    def fastafile                  = fasta                                                      ? "fastaFile = $\{PWD}${fasta}"                                                                 : ""
    def minimalcoverageperposition = task.ext.args?["general"]?["minimalcoverageperposition"]   ? "minimalCoveragePerPosition = ${task.ext.args["general"]["minimalcoverageperposition"]}"      : ""
    def minimalqualityperposition  = task.ext.args?["general"]?["minimalqualityperposition"]    ? "minimalQualityPerPosition = ${task.ext.args["general"]["minimalqualityperposition"]}"        : ""
    def shiftinquality             = task.ext.args?["general"]?["shiftinquality"]               ? "shiftInQuality = ${task.ext.args["general"]["shiftinquality"]}"                              : ""
    def snpfile                    = known_snps                                                 ? "SNPfile = \$PWD${known_snps}"                                                                : ""

    //"Target" configuration
    def target_bed                 = target_bed                                                 ? "captureRegions = ${target_bed}"                                                              : ""
    """
    touch config.txt

    echo "[general]" >> config.txt
    echo ${bedgraphoutput} >> config.txt
    echo ${breakpointthreshold} >> config.txt
    echo ${breakpointtype} >> config.txt
    echo "chrFiles = \${PWD}/${chr_directory.fileName}" >> config.txt
    echo "chrLenFile = \${PWD}/${chr_length.fileName}" >> config.txt
    echo ${coefficientofvariation} >> config.txt
    echo ${contamination} >> config.txt
    echo ${contaminationadjustment} >> config.txt
    echo ${degree} >> config.txt
    echo ${forcegccontentnormalization} >> config.txt
    echo ${gccontentprofile} >> config.txt
    echo ${mappability} >> config.txt
    echo ${intercept} >> config.txt
    echo ${mincnalength} >> config.txt
    echo ${minmappabilityperwindow} >> config.txt
    echo ${minexpectedgc} >> config.txt
    echo ${maxexpectedgc} >> config.txt
    echo ${minimalsubclonepresence} >> config.txt
    echo "maxThreads = ${task.cpus}" >> config.txt
    echo ${noisydata} >> config.txt
    echo "outputDir = ${prefix}" >> config.txt
    echo ${ploidy} >> config.txt
    echo ${printNA} >> config.txt
    echo ${readcountthreshold} >> config.txt
    echo ${sex} >> config.txt
    echo ${step} >> config.txt
    echo ${telocentromeric} >> config.txt
    echo ${uniquematch} >> config.txt
    echo ${window} >> config.txt

    echo "[sample]" >> config.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
