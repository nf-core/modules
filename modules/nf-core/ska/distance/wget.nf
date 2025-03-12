process WGET {
    input:
    val(url)
    val(outname)

    output:
    path("${outname}"), emit: outfile

    script:
    """
    wget -O ${outname} ${url}
    """
}
