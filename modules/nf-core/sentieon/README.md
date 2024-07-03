# How to use the sentieon modules

Except the sentieon bwaindex module, these nextflow sentieon modules require a license. Sentieon supply license in the form of a string-value (a url) or a file. It should be base64-encoded and stored in a nextflow secret named `SENTIEON_LICENSE_BASE64`. If a license string (url) is supplied, then the nextflow secret should be set like this:

```bash
nextflow secrets set SENTIEON_LICENSE_BASE64 $(echo -n <sentieon_license_string> | base64 -w 0)
```

If a license file is supplied, then the nextflow secret should be set like this:

```bash
nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
```

## Local Testing

``` bash
export SENTIEON_AUTH_MECH="GitHub Actions - token"
export SENTIEON_LICSRVR_IP=$(op read "op://Dev/Sentieon License Server/SENTIEON_LICSRVR_IP")
SENTIEON_ENCRYPTION_KEY=$(op read "op://Dev/Sentieon License Server/GitHub Secrets/SENTIEON_ENCRYPTION_KEY")
SENTIEON_LICENSE_MESSAGE=$(op read "op://Dev/Sentieon License Server/GitHub Secrets/SENTIEON_LICENSE_MESSAGE")
nextflow secrets set SENTIEON_AUTH_DATA $(python3 tests/modules/nf-core/sentieon/license_message.py encrypt --key "$SENTIEON_ENCRYPTION_KEY" --message "$SENTIEON_LICENSE_MESSAGE")
```

> [!NOTE]
> If this fails run `op signin` to flip to nf-core account
