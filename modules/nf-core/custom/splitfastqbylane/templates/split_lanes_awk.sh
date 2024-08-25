#!/bin/bash

if [[ "!{read1}" == *gz ]] ; then
    cat_="zcat"
else
    cat_="cat"
fi

function a() {
    awk \
        -v prefix=!{prefix} \
        -v readnumber=$1 \
        '
        BEGIN {FS = ":"}
        {
            lane=$(NF-3)
            flowcell=$(NF-4)
            outfastq=prefix"@"flowcell"_L00"lane"_R"readnumber".split.fastq.gz"
            print | "gzip > "outfastq
            for (i = 1; i <= 3; i++) {
                getline
                print | "gzip > "outfastq
            }
        }
        ' <( eval "$cat_ $2")
}

a 1 !{read1}
if [ ! -z !{read2} ] ; then
    a 2 !{read2}
fi

cat <<-END_VERSIONS > versions.yml
"!{task.process}":
    gawk: $(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
END_VERSIONS
