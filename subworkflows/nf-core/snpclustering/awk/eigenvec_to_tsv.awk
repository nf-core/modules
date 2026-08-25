NR==1 {
    sub(/^#/, "")
    if ($1 == "IID") {
        $1 = "sample_id"
    }
    print
    next
}
{ print }
