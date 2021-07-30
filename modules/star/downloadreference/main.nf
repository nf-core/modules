// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAR_DOWNLOADREFERENCE {
    tag 'genome_build'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'resources', meta:[:], publish_by_meta:[]) }

    input:
    val genome_build

    output:
    path "star-fusion-genome"       , emit: reference
    path "*.version.txt"            , emit: version

    script:
    def software        = getSoftwareName(task.process)

    if ( genome_build == "GRCh37" ){

    """
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz .
    mkdir -p star-fusion-genome
    tar -zxf GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz --strip-components=2 -C star-fusion-genome/
    rm GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz

    echo 'GRCh37_gencode_v19_CTAT_lib_Mar012021' > ${software}.version.txt
    """

    }
    else if (genome_build == "GRCh38"){

    """
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz .
    mkdir -p star-fusion-genome
    tar -zxf GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz --strip-components=2 -C star-fusion-genome/
    rm GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz

    echo 'GRCh38_gencode_v37_CTAT_lib_Mar012021' > ${software}.version.txt
    """

    }
    else{

    """
    echo \$genome_build > ${software}.version.txt
    mkdir -p star-fusion-genome
    cp ${software}.version.txt star-fusion-genome/
    """

    }
}
