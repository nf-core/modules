process SPOOF_SPIKES {

    input:
    path expression_matrix

    output:
    path 'spikes.txt'

    script:
    """
    head -n 50 $expression_matrix | \
        tail -n +2 | \
        awk '{print \$1}' > spikes.txt.tmp
    mv spikes.txt.tmp spikes.txt
    """
}
