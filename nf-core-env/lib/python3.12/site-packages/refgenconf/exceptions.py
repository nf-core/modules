""" Package exception types """

import abc

__all__ = [
    "DownloadJsonError",
    "GenomeConfigFormatError",
    "MissingAssetError",
    "MissingRecipeError",
    "MissingConfigDataError",
    "MissingGenomeError",
    "MissingSeekKeyError",
    "MissingTagError",
    "RefgenconfError",
    "UnboundEnvironmentVariablesError",
    "ConfigNotCompliantError",
    "RemoteDigestMismatchError",
]

DOC_URL = "http://refgenie.databio.org/en/latest/genome_config/"


class RefgenconfError(Exception):
    """Base exception type for this package"""

    __metaclass__ = abc.ABCMeta


class DownloadJsonError(RefgenconfError):
    """Non-OK response from a JSON download attempt"""

    def __init__(self, resp):
        super(DownloadJsonError, self).__init__(
            "No response provided" if resp is None else "JSON: {}".format(resp.json())
        )
        self.response = resp


class GenomeConfigFormatError(RefgenconfError):
    """Exception for invalid genome config file format."""

    def __init__(self, msg):
        spacing = " " if msg[-1] in ["?", ".", "\n"] else "; "
        suggest = "For config format documentation please see " + DOC_URL
        super(GenomeConfigFormatError, self).__init__(msg + spacing + suggest)


class MissingAssetError(RefgenconfError):
    """Error type for request of an unavailable genome asset."""

    pass


class MissingTagError(RefgenconfError):
    """Error type for request of an unavailable asset tag."""

    pass


class MissingSeekKeyError(RefgenconfError):
    """Error type for request of an unavailable asset seek key."""

    pass


class MissingRecipeError(RefgenconfError):
    """Error type for request of an unavailable recipe."""

    pass


class MissingConfigDataError(RefgenconfError):
    """Missing required configuration instance items"""

    pass


class ConfigNotCompliantError(GenomeConfigFormatError):
    """The format of the config file does not match required version/standards"""

    pass


class MissingGenomeError(RefgenconfError):
    """Error type for request of unknown genome/assembly."""

    pass


class UnboundEnvironmentVariablesError(RefgenconfError):
    """Use of environment variable that isn't bound to a value."""

    pass


class RemoteDigestMismatchError(RefgenconfError):
    """Remote digest of the parent asset does not match its local counterpart"""

    def __init__(self, asset, local_digest, remote_digest):
        msg = (
            "This asset is built from parent asset '{}', but for this parent, the remote does not "
            "match your local asset (local: {}; remote: {}). Refgenie will not pull this asset "
            "because the remote version was not built from the same parent asset you have locally.".format(
                asset, local_digest, remote_digest
            )
        )
        super(RemoteDigestMismatchError, self).__init__(msg)
