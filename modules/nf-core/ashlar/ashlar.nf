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
    params.output_channels = ''

    input:
        path file_in
    output:
        path '*.tif'

    script:
    println file_in
    def flip_x = params.flip_x == '' ? '' : '--flip-x'
    def flip_y = params.flip_y == '' ? '' : '--flip-y'
    def flip_mosaic_x = params.flip_mosaic_x == '' ? '' : '--flip-mosaic-x'
    def flip_mosaic_y = params.flip_mosaic_y == '' ? '' : '--flip-mosaic-y'
    def filter_sigma = params.filter_sigma == '' ? '' : '--filter-sigma ' + params.filter_sigma
    def plates = params.plates == '' ? '' : '--plates'
    def quiet = params.quiet == '' ? '' : '--quiet'
    def help = params.help == '' ? false : true
    def version = params.version == '' ? false : true
    def output_channels = params.output_channels == '' ? '' : '--output-channels ' + params.output_channels

    if(!help && !version) {
        """
        ashlar $file_in -o $params.output -c $params.channel $flip_x $flip_y $flip_mosaic_x $flip_mosaic_y -m $params.maximum_shift $filter_sigma --tile-size $params.tile_size $plates $quiet $output_channels
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
  ashlar(params.in) | wc | view { it.trim() }
}
