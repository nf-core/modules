process CONTROLFREEC_SOMATIC {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::control-freec=11.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/control-freec:11.6--h1b792b2_1':
        'quay.io/biocontainers/control-freec:11.6--h1b792b2_1' }"

    input:
    tuple val(meta), path(mpileup_normal), path(mpileup_tumor), path(cpn_normal), path(cpn_tumor), path(minipileup_normal), path(minipileup_tumor)
    path fasta
    path fai
    path snp_position
    path known_snps
    path known_snps_tbi
    path chr_length
    path chr_directory
    path mappability
    path target_bed

    output:
    tuple val(meta), path("*.txt"), emit: bedgraph, optional: true

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    //"General" configurations
    def bedgraphoutput              = task.ext.args?["general"]?["bedgraphoutput"]              ? "BedGraphOutput = ${task.ext.args["general"]["bedgraphoutput"]}"                              : ""
    //bedtools: not needed since pileup files are defined as input & we would need a new container, same for samtools and sambamba usage
    def chr_files                   = chr_directory                                             ? "chrFiles = \${PWD}/${chr_directory.fileName}"                                                : ""
    def chr_length                  = chr_length                                                ? "chrLength = \${PWD}/${chr_length.fileName}"                                                  : ""
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
    def matefile_normal             = mpileup_normal                                            ? "mateFile = \${PWD}/${mpileup_normal}"                                                        : ""
    def matecopynumberfile_normal   = cpn_normal                                                ? "mateCopyNumberFile = \${PWD}/${cpn_normal}"                                                  : ""
    def minipileup_normal           = minipileup_normal                                         ? "miniPileup = \${PWD}/${minipileup_normal}"                                                   : ""
    def inputformat_normal          = task.ext.args?["control"]?["inputformat"]                 ? "inputFormat = ${task.ext.args["control"]["inputformat"]}"                                    : ""
    def mateorientation_normal      = task.ext.args?["control"]?["mateorientation"]             ? "mateOrientation = ${task.ext.args["control"]["mateorientation"]}"                            : ""

    //"Sample" configuration
    def matefile_tumor             = mpileup_tumor                                              ? "mateFile = \${PWD}/${mpileup_tumor}"                                                         : ""
    def matecopynumberfile_tumor   = cpn_tumor                                                  ? "mateCopyNumberFile = \${PWD}/${cpn_tumor}"                                                   : ""
    def minipileup_tumor           = minipileup_tumor                                           ? "miniPileup = \${PWD}/${minipileup_tumor}"                                                    : ""
    def inputformat_tumor          = task.ext.args?["sample"]?["inputformat"]                   ? "inputFormat = ${task.ext.args["sample"]["inputformat"]}"                                     : ""
    def mateorientation_tumor      = task.ext.args?["sample"]?["mateorientation"]               ? "mateOrientation = ${task.ext.args["sample"]["mateorientation"]}"                             : ""

    //"BAF" configuration
    def makepileup                 = snp_position                                               ? "makePileup = \${PWD}/${snp_position}"                                                        : ""
    def fastafile                  = fasta                                                      ? "fastaFile = \${PWD}/${fasta}"                                                                : ""
    def minimalcoverageperposition = task.ext.args?["BAF"]?["minimalcoverageperposition"]       ? "minimalCoveragePerPosition = ${task.ext.args["BAF"]["minimalcoverageperposition"]}"          : ""
    def minimalqualityperposition  = task.ext.args?["BAF"]?["minimalqualityperposition"]        ? "minimalQualityPerPosition = ${task.ext.args["BAF"]["minimalqualityperposition"]}"            : ""
    def shiftinquality             = task.ext.args?["BAF"]?["shiftinquality"]                   ? "shiftInQuality = ${task.ext.args["BAF"]["shiftinquality"]}"                                  : ""
    def snpfile                    = known_snps                                                 ? "SNPfile = \$PWD/${known_snps}"                                                               : ""

    //"Target" configuration
    def target_bed                 = target_bed                                                 ? "captureRegions = ${target_bed}"                                                              : ""
    // TODO; fix output
    """
    touch config.txt

    echo "[general]" >> config.txt
    echo ${bedgraphoutput} >> config.txt
    echo ${breakpointthreshold} >> config.txt
    echo ${breakpointtype} >> config.txt
    echo ${chr_files} >> config.txt
    echo ${chr_length} >> config.txt
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
    echo "outputDir = ./" >> config.txt
    echo ${ploidy} >> config.txt
    echo ${printNA} >> config.txt
    echo ${readcountthreshold} >> config.txt
    echo ${sex} >> config.txt
    echo ${step} >> config.txt
    echo ${telocentromeric} >> config.txt
    echo ${uniquematch} >> config.txt
    echo ${window} >> config.txt

    echo "[control]" >> config.txt
    echo ${matefile_normal} >> config.txt
    echo ${matecopynumberfile_normal} >> config.txt
    echo ${minipileup_normal} >> config.txt
    echo ${inputformat_normal} >> config.txt
    echo ${mateorientation_normal} >> config.txt

    echo "[sample]" >> config.txt
    echo ${matefile_tumor} >> config.txt
    echo ${matecopynumberfile_tumor} >> config.txt
    echo ${minipileup_tumor} >> config.txt
    echo ${inputformat_tumor} >> config.txt
    echo ${mateorientation_tumor} >> config.txt

    echo "[BAF]" >> config.txt
    echo ${makepileup} >> config.txt
    echo ${fasta} >> config.txt
    echo ${minimalcoverageperposition} >> config.txt
    echo ${minimalqualityperposition} >> config.txt
    echo ${shiftinquality} >> config.txt
    echo ${snpfile} >> config.txt

    echo "[target]" >> config.txt
    echo ${target_bed} >> config.txt

    freec -conf config.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: \$(echo \$(freec -version 2>&1) | sed 's/^.*Control-FREEC  //; s/:.*\$//' | sed -e "s/Control-FREEC v//g" )
    END_VERSIONS
    """
}
