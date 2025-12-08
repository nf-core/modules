#!/usr/bin/env python
## Copyright 2018-2020 Stian Soiland-Reyes, The University of Manchester, UK
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.

"""
Parse arcp URIs.

Use is_arcp_uri() to detect of an URI string is using the 
arcp: URI scheme, in which case parse_arcp() can be used
to split it into its components.

The urlparse() function can be used as a replacement for
urllib.parse.urlparse() - supporting any URIs. If the URI is 
using the arcp: URI scheme, additional components are available
as from parse_arcp().
"""
__author__      = "Stian Soiland-Reyes <https://orcid.org/0000-0001-9842-9718>"
__copyright__   = "Copyright 2018-2020 The University of Manchester"
__license__     = "Apache License, version 2.0 (https://www.apache.org/licenses/LICENSE-2.0)"

from uuid import UUID, NAMESPACE_URL

try:
    import urllib.parse as urlp
except:
    import urlparse as urlp

from base64 import urlsafe_b64decode
from binascii import hexlify
import re

SCHEME="arcp"

def is_arcp_uri(uri):
    """Return True if the uri string uses the arcp scheme, otherwise False.
    """

    # tip: urllib will do lowercase for us
    return urlp.urlparse(uri).scheme == SCHEME

def parse_arcp(uri):
    """Parse an arcp URI string into its constituent parts.

    The returned object is similar to ``urllib.parse.urlparse()``
    in that it is a tuple of 
    ``(scheme,netloc,path,params,query,fragment)``
    with equally named properties, but it also adds
    properties for arcp fields:
    
    - prefix -- arcp authority prefix, e.g. "uuid", "ni" or "name", or None if prefix is missing
    - name -- arcp authority without prefix, e.g. "a4889890-a50a-4f14-b4e7-5fd83683a2b5" or "example.com"
    - uuid -- a ``uuid.UUID`` object if prefix is "uuid", otherwise None
    - ni -- the arcp alg-val value according to RFC6920 if prefix is "ni", otherwise None
    - hash -- the hash method and hash as a hexstring if prefix is "ni", otherwise None
    """

    return ARCPParseResult(*urlp.urlparse(uri))

def urlparse(uri):
    """Parse any URI string into constituent parts.

    The returned object is similar to 
    :func:`urllib.parse.urlparse()`
    in that it is a tuple of 
    ``(scheme,netloc,path,params,query,fragment)``
    with equally named properties, but if the 
    URI scheme is "arcp" this also adds
    arcp properties as in :func:`parse_arcp()`.
    """
    u = urlp.urlparse(uri)    
    if (u.scheme == SCHEME):
        return ARCPParseResult(*u)
    else:
        return u

class ARCPParseResult(urlp.ParseResult):
    """Result of parsing an arcp URI.

    This class does not detect if the arcp URI was valid
    according to the specification.

    This class extends :class:`urlllib.parse.ParseResult`
    adding arcp properties, some of which may be `None`.
    """
    __slots__ = ()

    def __init__(self, *args):
        if self.scheme != SCHEME:
            raise Exception("uri has scheme %s, expected %s" % 
                            (self.scheme, SCHEME))

    def _host_split(self):
        """Return (prefix,name) if authority has "," - 
        otherwise (None, authority).
        """
        if self.netloc and "," in self.netloc:
            return self.netloc.split(",", 1)
        else:
            return (None, self.netloc)

    @property
    def prefix(self):
        """The arcp prefix, e.g. "uuid", "ni", "name" or None if no prefix was present.
        """
        (prefix,name) = self._host_split()
        return prefix

    @property
    def name(self):
        """The URI's authority without arcp prefix.
        """
        (prefix,name) = self._host_split()
        return name
    
    @property
    def uuid(self):
        """The arcp UUID if the prefix is "uuid", otherwise None."""
        if self.prefix != "uuid":
            return None
        return UUID(self.name)
    
    @property
    def ni(self):
        """The arcp ni string if the prefix is "ni", otherwise None."""
        if self.prefix != "ni":
            return None
        if not _ALG_VAL.match(self.name):
            raise Exception("Invalid alg-val for ni, prefix: %s" % self.netloc)
        return self.name
    
    def ni_uri(self, authority=""):
        """The ni URI (RFC6920_) if the prefix is "ni", otherwise None.
        
        If the ``authority`` parameter is provided, 
        it will be used in the returned URI.

        .. _RFC6920: https://tools.ietf.org/search/rfc6920
        """
        ni = self.ni
        if ni is None:
            return None
        s = ("ni", authority, ni, None, None)
        return urlp.urlunsplit(s)


    def nih_uri(self):
        """The nih URI (RFC6920_) if the prefix is "ni", otherwise None.
        
        .. _RFC6920: https://tools.ietf.org/search/rfc6920
        """
        h = self.hash
        if h is None:
            return None
        (hash_method, hash_hex) = h
        
        segmented = _nih_segmented(hash_hex)
        checkdigit = _nih_checkdigit(hash_hex)

        path = "%s;%s;%s" % (hash_method, segmented, checkdigit)
        s = ("nih", None, path, None, None)
        return urlp.urlunsplit(s)
    
    def ni_well_known(self, base=""):
        """The ni .well-known URI (RFC5785_) if the prefix is 
        "ni", otherwise None.

        The parameter ``base``, if provided, should be an absolute URI like 
        ``"http://example.com/"`` - a relative URI is returned otherwise.

        .. _RFC5785: https://tools.ietf.org/html/rfc5785
        """
        (method, hash_b64) = self._ni_split()        
        if method is None:
            return None
        
        # .well-known is always at / (RFC5785)
        path = "/.well-known/ni/%s/%s" % (method, hash_b64)
        return urlp.urljoin(base, path)

    def _ni_split(self):
        """Split self.ni:
        """
        ni = self.ni
        if ni is None:
            return (None,None)
        # Already checked by self.ni regex
        #if not ";" in ni:
        #    raise Exception("invalid ni hash: %s" % ni)
        (method, hash_b64) = ni.split(";", 1)
        return (method, hash_b64)
        

    @property
    def hash(self):
        """A tuple (hash_method,hash_hex) if the prefix is "ni", 
        otherwise None.
        """
        (method, hash_b64) = self._ni_split()
        if method is None:
            return None
        # re-instate padding as urlsafe_base64decode is strict
        missing_padding = 4 - (len(hash_b64) % 4)
        hash_b64 += "=" * missing_padding
        hash_bytes = urlsafe_b64decode(hash_b64)
        hash_hex = hexlify(hash_bytes).decode("ascii")
        return (method.lower(), hash_hex)
    
    def __repr__(self):
        props = ["scheme='arcp'"]
        props += ["prefix='%s'" % self.prefix or ""]
        props += ["name='%s'" % self.name or ""]

        if self.uuid is not None:
            props += ["uuid=%s" % self.uuid]
        if self.ni is not None:
            props += ["ni='%s'" % self.ni]
            # Avoid Exception in __repr__
            if ";" in self.ni:
                props += ["hash=('%s', '%s'" % self.hash]

        # Traditional URI properties
        props += ["path='%s'" % self.path or ""]
        props += ["query='%s'" % self.query or ""]
        props += ["fragment='%s'" % self.fragment or ""]
        return "ARCPParseResult(%s)" % ",".join(props)

    def __str__(self):
        return self.geturl()

def _alg_val_regex():
    """Compile regular expression for RFC6920_ alg-val production

    .. _RFC6920: https://www.ietf.org/rfc/rfc6920
    """
    # unreserved    = ALPHA / DIGIT / "-" / "." / "_" / "~"
    unreserved = r"[A-Za-z0-9-._~]"
    # alg            = 1*unreserved
    alg = r"(" + unreserved + r"+)"
    # val            = 1*unreserved
    val = r"(" + unreserved + r"+)"
    # alg-val        = alg ";" val
    alg_val = r"^" + alg + ";" + val + r"$"
    return re.compile(alg_val)
_ALG_VAL = _alg_val_regex()

def _nih_segmented(h, grouping=6):
    """Segment hex-hash with dashes in nih style RFC6920_

    >>> _nih_segmented("0123456789abcdef")
    "012345-6789ab-cdef"

    .. _RFC6920: https://www.ietf.org/rfc/rfc6920
    """
    segmented = []
    while h:
        segmented.append(h[:grouping])
        h = h[grouping:]
    return "-".join(segmented)

def _nih_checkdigit(h):
    """Luhn mod N algorithm in base 16 (hex) according to RFC6920_

    .. _RFC6920: https://www.ietf.org/rfc/rfc6920
    """
    ## Adopted from https://en.wikipedia.org/wiki/Luhn_mod_N_algorithm 
    ## pseudocode
    factor = 2
    total = 0
    base = 16
    digits = len(h)        
    # 0 if digits has even length, 1 if odd
    # (as we start doubling with the very last digit)
    parity = digits % 2
    for x in range(digits):
        digit = int(h[x], 16)            
        if x % 2 != parity:
            # double every second digit
            digit *= 2
            # slight less efficient, but more verbose:
            # if > 16:
            #   total += digit - 16 + 1
            # else: 
            #   total + digit
            total += sum(divmod(digit, 16))
        else:
            # Not doubled, must be <16
            total += digit
    # checkdigit that needs to be added to total
    # to get 0 after modulus
    remainder = (16-total) % 16
    # Return as hex digit
    return "%x" % remainder
