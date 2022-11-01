params.str = 'Hello world!'

process splitLetters {
  output:
    path 'chunk_*'

  """
  printf '${params.str}' | split -b 6 - chunk_
  """
}

process convertToUpper {
  input:
    path x
  output:
    stdout

  """
  cat $x | tr '[a-z]' '[A-Z]'
  """
}

process ashlar_works_1 {
    conda 'bioconda::ashlar=1.17.0'

    input:
        path file_in
    output:
        stdout

    """
    ashlar $file_in -o foo.ome.tif
    """

    /*
    output:
        stdout
    """
    ashlar --help
    """
    */
}

process ashlar_works_2 {
    conda 'bioconda::ashlar=1.17.0'

    input:
        path file_in
    output:
        stdout

    """
    ashlar $file_in -o $params.out
    """

}

process wc {
    input:
        path file_in
    output:
        stdout

    """
        wc -l $file_in
    """
}

process ashlar {
    conda 'bioconda::ashlar=1.17.0'
    params.out = 'ashlar_output.ome.tif'

    input:
        path file_in
    output:
        path '*.ome.tif'

    """
    ashlar $file_in -o $params.out
    """

}

workflow {
  // splitLetters | flatten | convertToUpper | view { it.trim() }
  // ashlar(params.in) | view { it.trim() }
  //wc(params.in) | view { it.trim() }
  ashlar(params.in) | wc | view { it.trim() }
}
