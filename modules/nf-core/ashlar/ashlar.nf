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
    params.output = 'ashlar_output.ome.tif'
    params.channel = '0'
    params.flip_x = ''
    params.flip_y = ''
    params.flip_mosaic_x = ''
    params.flip_mosaic_y = ''
    params.maximum_shift = '15'
    params.filter_sigma = ''
    params.tile_size = '1024'
    params.plates = ''
    params.quiet = ''
    params.help = ''
    params.version = ''

    input:
        path file_in
    output:
        path '*.tif'

    script:
    def flip_x = params.flip_x == '' ? '' : '--flip-x'
    def flip_y = params.flip_y == '' ? '' : '--flip-y'
    def flip_mosaic_x = params.flip_mosaic_x == '' ? '' : '--flip-mosaic-x'
    def flip_mosaic_y = params.flip_mosaic_y == '' ? '' : '--flip-mosaic-y'
    def filter_sigma = params.filter_sigma == '' ? '' : '--filter-sigma ' + params.filter_sigma
    def plates = params.plates == '' ? '' : '--plates'
    def quiet = params.quiet == '' ? '' : '--quiet'
    def help = params.help == '' ? false : true
    def version = params.version == '' ? false : true

    if(!help && !version) {
        """
        ashlar $file_in -o $params.output -c $params.channel $flip_x $flip_y $flip_mosaic_x $flip_mosaic_y -m $params.maximum_shift $filter_sigma --tile-size $params.tile_size $plates $quiet
        """
    }
    else if(help && version) {
        """
        echo 'Please select either the help or version input parameter, but not both.'
        """
    }
    else if(help) {
        """
        ashlar --help
        """
    }
    else if(version) {
        """
        ashlar --version
        """
    }
    else {
        """
        ashlar --help
        """
    }
}

workflow {
  // splitLetters | flatten | convertToUpper | view { it.trim() }
  // ashlar(params.in) | view { it.trim() }
  //wc(params.in) | view { it.trim() }
  ashlar(params.in) | wc | view { it.trim() }
}
