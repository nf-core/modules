

process MTMALIGN_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5bcf71dc66dac33d8e003c5e78043b80f5c7f269:8f0e486d46f7ab38892c1a8f78d2894a549d03b5-0':
        'biocontainers/mulled-v2-5bcf71dc66dac33d8e003c5e78043b80f5c7f269:8f0e486d46f7ab38892c1a8f78d2894a549d03b5-0' }"

    input:
    tuple val(meta), path(pdbs)
    val(compress)

    output:
    tuple val(meta), path("${prefix}.aln${compress ? '.gz' : ''}"), emit: alignment
    tuple val(meta), path("${prefix}.pdb${compress ? '.gz' : ''}"), emit: structure
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // mTMalign is not capable of writing to stdout
    // if -o /dev/stdout is specified, the output file will be polluted with debug messages emitted by mTMalign
    """
    # decompress input files if required
    if ls ./*.pdb.gz 2&> /dev/null; then # check if any files are compressed; calling unpigz with an empty arg will cause it to panic
        unpigz -d ./*.pdb.gz
    fi

    # construct input file for mtmalign
    ls *.pdb | sed s/\\ /\\n/ > input_list.txt

    mtm-align -i input_list.txt -o ${prefix}.pdb
    # -o does not affect the fasta naming, so move it to the new name
    mv ./mTM_result/result.fasta ${prefix}.aln
    mv ./mTM_result/${prefix}.pdb ${prefix}.pdb
    # Remove ".pdb" from the ids in the alignment file
    sed -i 's/\\.pdb//g' ${prefix}.aln

    # compress both output files
    if ${compress}; then
        pigz -p ${task.cpus} ${prefix}.aln ${prefix}.pdb
    fi

    # mtm-align -v prints the wrong version 20180725, so extract it from the cosmetic output in the help message
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mTM-align: \$( mtm-align -h | grep -e "\\(Version [[:digit:]]*\\)" | grep -oe "[[:digit:]]*" )
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln${compress ? '.gz' : ''}
    touch ${prefix}.pdb${compress ? '.gz' : ''}

    # mtm-align -v prints the wrong version 20180725, so extract it from the cosmetic output in the help message
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mTM-align: \$( mtm-align -h | grep -e "\\(Version [[:digit:]]*\\)" | grep -oe "[[:digit:]]*" )
        pigz: \$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' ))
    END_VERSIONS
    """
}
