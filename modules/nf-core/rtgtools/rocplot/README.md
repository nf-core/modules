# RTG Tools - rocplot

[RTG Tools GitHub link](https://github.com/RealTimeGenomics/rtg-tools)

## Overview

Plot ROC curves from vcfeval ROC data files, either to an image, or an interactive GUI.

## Command help

`rtg rocplot --help`

```bash

Usage: rtg rocplot [OPTION]... FILE+
                   [OPTION]... --curve STRING

Plot ROC curves from vcfeval ROC data files, either to an image, or an interactive GUI.

File Input/Output
      --curve=STRING          ROC data file with title optionally specified (path[=title]). May be specified 0 or more times
      --png=FILE              if set, output a PNG image to the given file
      --svg=FILE              if set, output a SVG image to the given file
      --zoom=STRING           show a zoomed view with the given coordinates, supplied in the form <xmax>,<ymax> or <xmin>,<ymin>,<xmax>,<ymax>
      FILE+                   ROC data file. May be specified 0 or more times

Reporting
      --cmd=FILE              if set, print rocplot command used in previously saved image
      --hide-sidepane         if set, hide the side pane from the GUI on startup
      --interpolate           if set, interpolate curves at regular intervals
      --line-width=INT        sets the plot line width (Default is 2)
      --palette=STRING        name of color palette to use. Allowed values are [blind_13, blind_15, blind_8, brewer_accent, brewer_dark2, brewer_paired,
                              brewer_pastel1, brewer_pastel2, brewer_set1, brewer_set2, brewer_set3, classic] (Default is classic)
      --plain                 if set, use a plain plot style
  -P, --precision-sensitivity if set, plot precision vs sensitivity rather than ROC
      --scores                if set, show scores on the plot
  -t, --title=STRING          title for the plot

Utility
  -h, --help                  print help on command-line flag usage
```
