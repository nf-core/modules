nextflow.enable.dsl=2

process ZERO_UUID {

    input:
    val(file_in)
    val(offset)

    when:
    file_in != "versions.yml"

    script:
    """

    IFS=',' read -ra FILE_IN_ARRAY <<< "$file_in"
    FILE_PATH_REMOVED_PREFIX="/\${FILE_IN_ARRAY[2]#*/}"
    FILE_PATH=\${FILE_PATH_REMOVED_PREFIX::-1}

    echo -n "00000000-0000-0000-0000-000000000000" | dd of=\${FILE_PATH} bs=1 seek=$offset conv=notrunc

    """

}
