nextflow.enable.dsl=2

process ZERO_UUID {

    input:
    val(file_in)
    val(offset)

    when:
    file_in != "versions.yml"

    script:
    def file_path = file_in[1]

    """

    echo -n "00000000-0000-0000-0000-000000000000" | dd of=$file_path bs=1 seek=$offset conv=notrunc

    """

}
