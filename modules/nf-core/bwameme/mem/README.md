# bwameme/mem options

bwa-meme is a high-throughput aligner that produces alignments faster than samtools can consume them. To prevent bwa-meme from stalling while samtools processes data, mbuffer is used as an intermediary to absorb alignments until samtools is ready. The mbuffer and samtools options are configurable via `ext.args2` and `ext.args3` respectively.

## Configuring mbuffer (`ext.args2`)

`ext.args2` is passed to mbuffer. The default buffer size is 3GB (`-m 3072M`). If `-m` is not present in `ext.args2`, the default is injected automatically.

The mbuffer size should match the total memory allocated to `samtools sort` (`-m` × `-@`) so it can absorb bwa-meme output while samtools is flushing its sort buffer to disk.

```
withName: 'BWAMEME_MEM' {
    ext.args2 = '-m 20480M'  // 20GB mbuffer to match samtools sort total (e.g. -m 1024M x -@ 20)
}
```

## Configuring samtools (`ext.args3`)

`ext.args3` is passed to `samtools sort` or `samtools view` depending on the `sort_bam` input. Defaults are injected if not supplied:

- `-@ 3` (threads) — always injected if `-@` is absent
- `-m 1024M` (memory per thread) — injected if `-m` is absent and `sort_bam` is true

```
withName: 'BWAMEME_MEM' {
    ext.args3 = '-@ 20 -m 1024M'  // 20 threads, 1GB per thread = 20GB total
}
```

## Example: tuning for human genome alignment

For a large reference (e.g. human genome) on a well-resourced machine, you may want to increase both values. The mbuffer size should equal the total samtools sort memory:

```
withName: 'BWAMEME_MEM' {
    ext.args2 = '-m 20480M'   // mbuffer = samtools total (20 threads x 1GB)
    ext.args3 = '-@ 20 -m 1024M'
}
```

## CRAM output

To produce CRAM output, pass `--output-fmt` via `ext.args3`:

```
withName: 'BWAMEME_MEM' {
    ext.args3 = '--output-fmt cram'
}
```

A FASTA reference must be provided as `input[2]` when using CRAM output.
