BEGIN {
    FS = OFS = "\t"
}

NR == 1 {
    sub(/^#/, "")
    $0 = $0
    for (i = 1; i <= NF; i++) {
        if ($i == "IID") {
            $i = "sample_id"
        }
    }
    print
    next
}

{
    print
}
