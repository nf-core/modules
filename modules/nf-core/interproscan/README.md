# Updating the docker container and making a new module release

As of March 2025, the latest version of [InterProScan on Bioconda](https://bioconda.github.io/recipes/interproscan/README.html) is 5.59_91, which is from [October 17, 2022](https://github.com/ebi-pf-team/interproscan/releases/tag/5.59-91.0). This Dockerfile builds a new version of InterProScan.

1. Edit the Dockerfile

Update the InterProScan versions in this line in the Dockerfile:

```bash
ENV INTERPROSCAN_VER=<VERSION>
```

2. Create and test the container:

Make sure to `export INTERPROSCAN_VER=<VERSION>` so that the build uses the right tags. If you're on arm64 architecture (e.g. Apple Silicon), you may need to run `export DOCKER_DEFAULT_PLATFORM=linux/amd64` to "force" the architecture to be `x86_64`/`amd64`.

You can do `make build` from the Makefile (requires `export INTERPROSCAN_VER=<VERSION>`) or:

```bash
docker build . -t quay.io/nf-core/interproscan:<VERSION>
```

3. Access rights are needed to push the container to the Dockerhub nfcore organization, please ask a core team member to do so.

You can do `make push` from the Makefile (requires `export INTERPROSCAN_VER=<VERSION>`) or:

```bash
docker push quay.io/nf-core/interproscan:<VERSION>
```

## Test Data Output

### `test.gff3.gz`

```
$ gunzip -c test.gff3.gz
##gff-version 3
##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/sofa.obo?revision=1.269
##interproscan-version 5.73-104.0
##sequence-region test_1 1 409
test_1  .       polypeptide     1       409     .       +       .       ID=test_1;md5=e622960f07c2857c8ac01636880665b0
test_1  Coils   protein_match   303     323     .       +       .       date=25-03-2025;Target=test_1 303 323;ID=match$1_303_323;Name=Coil;status=T
##FASTA
>test_1
MSSHPIQVFSEIGKLKKVMLHRPGKELENLLPDYLERLLFDDIPFLEDAQKEHDAFAQAL
RDEGIEVLYLEQLAAESLTSQEIRDQFIEEYLDEANIRDRHTKVAIRELLHGIKDNQELV
EKTMAGIQKVELPEIPDEAKDLTDLVESDYPFAIDPMPNLYFTRDPFATIGNAVSLNHMF
ADTRNRETLYGKYIFKYHPIYGGKVDLVYNREEDTRIEGGDELVLSKDVLAVGISQRTDA
ASIEKLLVNIFKKNVGFKKVLAFEFANNRKFMHLDTVFTMVDYDKFTIHPEIEGDLHVYS
VTYENEKLKIVEEKGDLAELLAQNLGVEKVHLIRCGGGNIVAAAREQWNDGSNTLTIAPG
VVVVYDRNTVTNKILEESGLRLIKIRGSELVRGRGGPRCMSMPFEREEV
>match$1_303_323
YENEKLKIVEEKGDLAELLAQ
```

### `test.json.gz`

```
$ gunzip -c test.json.gz
{
  "interproscan-version": "5.73-104.0",
  "results": [ {
    "sequence" : "MSSHPIQVFSEIGKLKKVMLHRPGKELENLLPDYLERLLFDDIPFLEDAQKEHDAFAQALRDEGIEVLYLEQLAAESLTSQEIRDQFIEEYLDEANIRDRHTKVAIRELLHGIKDNQELVEKTMAGIQKVELPEIPDEAKDLTDLVESDYPFAIDPMPNLYFTRDPFATIGNAVSLNHMFADTRNRETLYGKYIFKYHPIYGGKVDLVYNREEDTRIEGGDELVLSKDVLAVGISQRTDAASIEKLLVNIFKKNVGFKKVLAFEFANNRKFMHLDTVFTMVDYDKFTIHPEIEGDLHVYSVTYENEKLKIVEEKGDLAELLAQNLGVEKVHLIRCGGGNIVAAAREQWNDGSNTLTIAPGVVVVYDRNTVTNKILEESGLRLIKIRGSELVRGRGGPRCMSMPFEREEV",
    "md5" : "e622960f07c2857c8ac01636880665b0",
    "matches" : [ {
        "signature" : {
        "accession" : "Coil",
        "name" : "Coil",
        "description" : null,
        "signatureLibraryRelease" : {
            "library" : "COILS",
            "version" : "2.2.1"
        },
        "entry" : null
        },
        "locations" : [ {
        "start" : 303,
        "end" : 323,
        "representative" : false,
        "location-fragments" : [ {
            "start" : 303,
            "end" : 323,
            "dc-status" : "CONTINUOUS"
        } ]
        } ],
        "model-ac" : "Coil"
    } ],
    "xref" : [ {
        "name" : "test_1",
        "id" : "test_1"
    } ]
    } ]
}
```

### `test.tsv.gz`

```
$ gunzip -c test.tsv.gz
test_1  e622960f07c2857c8ac01636880665b0        409     Coils   Coil    Coil    303     323     -       T       25-03-2025      -       -       -       -
```

### `test.xml.gz`

```
$ gunzip -c test.xml.gz
<?xml version="1.0" encoding="UTF-8"?><protein-matches xmlns="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas" interproscan-version="5.73-104.0">
  <protein>
    <sequence md5="e622960f07c2857c8ac01636880665b0">MSSHPIQVFSEIGKLKKVMLHRPGKELENLLPDYLERLLFDDIPFLEDAQKEHDAFAQALRDEGIEVLYLEQLAAESLTSQEIRDQFIEEYLDEANIRDRHTKVAIRELLHGIKDNQELVEKTMAGIQKVELPEIPDEAKDLTDLVESDYPFAIDPMPNLYFTRDPFATIGNAVSLNHMFADTRNRETLYGKYIFKYHPIYGGKVDLVYNREEDTRIEGGDELVLSKDVLAVGISQRTDAASIEKLLVNIFKKNVGFKKVLAFEFANNRKFMHLDTVFTMVDYDKFTIHPEIEGDLHVYSVTYENEKLKIVEEKGDLAELLAQNLGVEKVHLIRCGGGNIVAAAREQWNDGSNTLTIAPGVVVVYDRNTVTNKILEESGLRLIKIRGSELVRGRGGPRCMSMPFEREEV</sequence>
    <xref id="test_1" name="test_1"/>
    <matches>
      <coils-match>
        <signature ac="Coil" name="Coil">
          <signature-library-release library="COILS" version="2.2.1"/>
        </signature>
        <model-ac>Coil</model-ac>
        <locations>
          <coils-location start="303" end="323" representative="false">
            <location-fragments>
              <coils-location-fragment start="303" end="323" dc-status="CONTINUOUS"/>
            </location-fragments>
          </coils-location>
        </locations>
      </coils-match>
    </matches>
  </protein>
</protein-matches>
```
