"""auth exceptions"""


class PEPHubResponseException(Exception):
    """Request response exception. Used when response != 200"""

    def __init__(self, reason: str = ""):
        """
        Optionally provide explanation for exceptional condition.
        :param str reason: some context or perhaps just a value that
            could not be interpreted as an accession
        """
        super(PEPHubResponseException, self).__init__(reason)


class PEPHubTokenExchangeException(Exception):
    """Exception in exchanging device code on token == 400"""

    def __init__(self, reason: str = ""):
        """
        Optionally provide explanation for exceptional condition.

        :param str reason: some context or perhaps just a value that
            could not be interpreted as an accession
        """
        super(PEPHubTokenExchangeException, self).__init__(reason)
