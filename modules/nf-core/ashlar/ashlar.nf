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
    println 'running process ashlar in ashlar.nf'
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
    params.ffp = ''
    params.dfp = ''

    input:
        val file_in
        val file_ffp
        val file_dfp
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
    def output_channels = params.output_channels == '' ? '' : '--output-channels ' + params.output_channels
    def ffp = params.ffp == '' ? '' : '--ffp ' + file_ffp
    def dfp = params.dfp == '' ? '' : '--dfp ' + file_dfp

    if(!help && !version) {
        """
        ashlar $file_in -o $params.output -c $params.channel $flip_x $flip_y $flip_mosaic_x $flip_mosaic_y -m $params.maximum_shift $filter_sigma --tile-size $params.tile_size $plates $quiet $output_channels $ffp $dfp
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
    input_file_list = params.in.split(',').join(' ')
    ffp_file_list = []
    if(params.ffp) {
        ffp_file_list = params.ffp.split(',').join(' ')
    }
    dfp_file_list = []
    if(params.dfp) {
        dfp_file_list = params.dfp.split(',').join(' ')
    }

    ashlar(input_file_list, ffp_file_list, dfp_file_list) | wc | view { it.trim() }
}
