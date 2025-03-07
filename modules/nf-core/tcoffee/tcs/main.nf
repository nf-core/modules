process TCOFFEE_TCS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/t-coffee_pigz:91ac7e26b23bb246':
        'community.wave.seqera.io/library/t-coffee_pigz:7d1373a24f76afe6' }"

    input:
    tuple val(meta) , path(msa)
    tuple val(meta2), path(lib)

    output:
    tuple val(meta), path("*.tcs")   , emit: tcs
    tuple val(meta), path("*.scores"), emit: scores
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def lib_arg       = lib ? "-lib ${lib}" : ""
    def header        = meta.keySet().join(",")
    def values        = meta.values().join(",")
    def unzipped_name = msa.toString() - '.gz'
    """
    export TEMP='./'
    export TMP_4_TCOFFEE="./"
    export HOME="./"

    filename=${msa}
    if [[ \$(basename $msa) == *.gz ]] ; then
        unpigz -f $msa
        filename=${unzipped_name}
    fi

    # Bad hack to circumvent t_coffee bug
    # Issue described already in: https://github.com/cbcrg/tcoffee/issues/3
    # Add an A in front of filename if the file begins with A
    first_letter_filename=\${filename:0:1}
    if [ "\$first_letter_filename" == "A" ]; then input="A"\$filename; else input=\$filename;  fi

    t_coffee -infile \$input \
        -evaluate -output=score_ascii \
        ${lib_arg} \
        -outfile ${prefix}.tcs

    # Add metadata info to output file
    echo "${header},TCS" > "${prefix}.scores"

    # Add values
    scores=\$(grep 'SCORE=' ${prefix}.tcs | cut -d '=' -f 2 )
    echo "${values},\$scores" >> "${prefix}.scores"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    export TMP_4_TCOFFEE="./"
    export HOME="./"

    touch ${prefix}.tcs
    touch ${prefix}.scores

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
