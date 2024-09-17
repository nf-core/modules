# Sentieon Modules

## Usecases we're trying to cover

In rank order:

1. Paying Sentieon users running a license server inside their AWS account(or on-prem).
2. GitHub actions using a _test_ license server that can receive traffic from the internet. This has an auth layer to keep anyone from connecting to it.
3. Pipeline developers working with Sentieon modules locally with their own license server.
4. Pipeline developers working with Sentieon modules locally with access to the `SENTIEON_LICSRVR_IP` and `SENTIEON_ENCRYPTION_KEY`.
5. People using a local file(Sentieon doesn't use these test licenses often anymore.)

## How to use the sentieon modules

Except the Sentieon bwaindex module, these nextflow Sentieon modules require a license. Sentieon supply license in the form of a string-value (a url ex: `12.12.12.12:8990`) or a file. The file should be base64-encoded and stored in a nextflow secret named `SENTIEON_LICENSE_BASE64`. If a license string (url) is supplied, then the nextflow secret should be set like this:

### Most Users

<!-- NOTE Might be SENTIEON_LICENSE = "$SENTIEON_LICSRVR_IP" -->

```nextflow title="nextflow.config"
env {
    SENTIEON_LICSRVR_IP = "$SENTIEON_LICSRVR_IP"
    // or if you really want to use `nextflow secrets`
    SENTIEON_LICSRVR_IP = secrets.SENTIEON_LICSRVR_IP
}
```

### License File Users

If a license file is supplied, then the nextflow secret should be set like this:

```bash
nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
```

And in the config set

```nextflow title="nextflow.config"
env {
    SENTIEON_LICENSE_BASE64 = secrets.SENTIEON_LICENSE_BASE64
}
```

## GitHub Actions

### Spiking in all the Sentieon variables

Previously we tried to get everything into one Nextflow secret. This made the setup extremely difficult to follow and reproduce.

The issue is that Nextflow secrets that are declared in a process cause it to fail if they're not used. This made it so we had to branch our logic off the size of the secret. Again, just making a confusing setup, even more confusing.

We tried to follow as closely as possible [DonFreed/docker-actions-test](https://github.com/DonFreed/docker-actions-test) and keep the process simple. At the end of the day Sentieon is also validating the license servers and can just make a new license for us if someone tries to misuse ours.

The main experience should be on end users of the Sentieon modules ease of use.

### Threat Models

After working with [@DonFreed](https://github.com/DonFreed) on a [rework for GitHub actions](https://github.com/nf-core/modules/pull/5856), we determined that using Nextflow secrets wasn't really necessary outside of the `SENTIEON_AUTH_DATA`.

Because the SENTIEON_ENCRYPTION_KEY and SENTIEON_LICENSE_MESSAGE are both stored in GitHub Secrets(and 1Password) contributors will never have access to them. They are just used to interact with the `license_message` script that generates some auth data to then come into contact with the server.

The `SENTIEON_AUTH_DATA` is only valid for 1 day, and the license that we have is valid for only a few CPUs.

The Server IP doesn't matter either because they would also need the `SENTIEON_ENCRYPTION_KEY` to authenticate with the license.

## Local Testing

```bash
export SENTIEON_LICSRVR_IP=$(op read "op://Dev/Sentieon License Server/SENTIEON_LICSRVR_IP")
```

For `SENTIEON_LICSRVR_IP` you can use your own server IP.

Additionally, the following may be necessary. If you don't know either key, you probably don't need it:

<details markdown="1">
<summary>Optional configuration</summary>

```bash
export SENTIEON_AUTH_MECH="GitHub Actions - token"
SENTIEON_ENCRYPTION_KEY=$(op read "op://Dev/Sentieon License Server/GitHub Secrets/SENTIEON_ENCRYPTION_KEY")
SENTIEON_LICENSE_MESSAGE=$(op read "op://Dev/Sentieon License Server/GitHub Secrets/SENTIEON_LICENSE_MESSAGE")
nextflow secrets set SENTIEON_AUTH_DATA $(python3 tests/modules/nf-core/sentieon/license_message.py encrypt --key "$SENTIEON_ENCRYPTION_KEY" --message "$SENTIEON_LICENSE_MESSAGE")
```

> [!NOTE]
> If this fails run `op signin` to flip to nf-core account

</details>
