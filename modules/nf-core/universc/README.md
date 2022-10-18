# UniverSC

## Single-cell processing across technologies

UniverSC is an open-source single-cell pipeline that runs across platforms on various technologies.

## Maintainers

Tom Kelly (RIKEN, IMS)

Kai Battenberg (RIKEN CSRS/IMS)

Contact: <first name>.<family name>[at]riken.jp

## Implementation

This container runs Cell Ranger v3.0.2 installed from source on MIT License on GitHub with
modifications for compatibility with updated dependencies. All software is installed from
open-source repositories and available for reuse.

It is _not_ subject to the 10X Genomics End User License Agreement (EULA).
This version allows running Cell Ranger v3.0.2 on data generated from any experimental platform
without restrictions. However, updating to newer versions on Cell Ranger subject to the
10X EULA is not possible without the agreement of 10X Genomics.

To comply with licensing and respect 10X Genomics Trademarks, the 10X Genomics logo
has been removed from HTML reports, the tool has been renamed, and proprietary
closed-source tools to build Cloupe files are disabled.

It is still suffient to generate summary reports and count matrices compatible with
single-cell analysis tools available for 10X Genomics and Cell Ranger output format
in Python and R packages.

## Disclaimer

We are third party developers not affiliated with 10X Genomics or any other vendor of
single-cell technologies. We are releasing this code on an open-source license which calls Cell Ranger
as an external dependency.

## Licensing

This package is provided open-source on a GPL-3 license. This means that you are free to use and
modify this code provided that they also contain this license.

## Updating the package

The tomkellygenetics/universc:<VERSION> container is automatically updated with tomkellygenetics/universc:latest.

To build an updated container use the Dockerfile provided here:

[https://github.com/minoda-lab/universc/blob/master/Dockerfile](https://github.com/minoda-lab/universc/blob/master/Dockerfile)
