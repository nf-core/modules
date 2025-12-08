from __future__ import annotations

import calendar
import contextlib
import hashlib
import io
import os
import platform
import shutil
import tempfile
import time
import typing as t

import requests

_LASTMOD_FMT = "%a, %d %b %Y %H:%M:%S %Z"


def _base_cache_dir() -> str | None:
    sysname = platform.system()

    # on windows, try to get the appdata env var
    # this *could* result in cache_dir=None, which is fine, just skip caching in
    # that case
    if sysname == "Windows":
        cache_dir = os.getenv("LOCALAPPDATA", os.getenv("APPDATA"))
    # macOS -> app support dir
    elif sysname == "Darwin":
        cache_dir = os.path.expanduser("~/Library/Caches")
    # default for unknown platforms, namely linux behavior
    # use XDG env var and default to ~/.cache/
    else:
        cache_dir = os.getenv("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))

    return cache_dir


def _resolve_cache_dir(dirname: str) -> str | None:
    cache_dir = _base_cache_dir()
    if cache_dir:
        cache_dir = os.path.join(cache_dir, "check_jsonschema", dirname)
    return cache_dir


def _lastmod_from_response(response: requests.Response) -> float:
    try:
        return calendar.timegm(
            time.strptime(response.headers["last-modified"], _LASTMOD_FMT)
        )
    # OverflowError: time outside of platform-specific bounds
    # ValueError: malformed/unparseable
    # LookupError: no such header
    except (OverflowError, ValueError, LookupError):
        return 0.0


def _get_request(
    file_url: str, *, response_ok: t.Callable[[requests.Response], bool]
) -> requests.Response:
    num_retries = 2
    r: requests.Response | None = None
    for _attempt in range(num_retries + 1):
        try:
            r = requests.get(file_url, stream=True)
        except requests.RequestException as e:
            if _attempt == num_retries:
                raise FailedDownloadError("encountered error during download") from e
            continue
        if r.ok and response_ok(r):
            return r
    assert r is not None
    raise FailedDownloadError(
        f"got response with status={r.status_code}, retries exhausted"
    )


def _atomic_write(dest: str, content: bytes) -> None:
    # download to a temp file and then move to the dest
    # this makes the download safe if run in parallel (parallel runs
    # won't create a new empty file for writing and cause failures)
    fp = tempfile.NamedTemporaryFile(mode="wb", delete=False)
    fp.write(content)
    fp.close()
    shutil.copy(fp.name, dest)
    os.remove(fp.name)


def _cache_hit(cachefile: str, response: requests.Response) -> bool:
    # no file? miss
    if not os.path.exists(cachefile):
        return False

    # compare mtime on any cached file against the remote last-modified time
    # it is considered a hit if the local file is at least as new as the remote file
    local_mtime = os.path.getmtime(cachefile)
    remote_mtime = _lastmod_from_response(response)
    return local_mtime >= remote_mtime


def url_to_cache_filename(ref_url: str) -> str:
    """
    Given a schema URL, convert it to a filename for caching in a cache dir.

    Rules are as follows:
    - the base filename is an sha256 hash of the URL
    - if the filename ends in an extension (.json, .yaml, etc) that extension
      is appended to the hash

    Preserving file extensions preserves the extension-based logic used for parsing, and
    it also helps a local editor (browsing the cache) identify filetypes.
    """
    filename = hashlib.sha256(ref_url.encode()).hexdigest()
    if "." in (last_part := ref_url.rpartition("/")[-1]):
        _, _, extension = last_part.rpartition(".")
        filename = f"{filename}.{extension}"
    return filename


class FailedDownloadError(Exception):
    pass


class CacheDownloader:
    def __init__(self, cache_dir: str, *, disable_cache: bool = False) -> None:
        self._cache_dir = _resolve_cache_dir(cache_dir)
        self._disable_cache = disable_cache

    def _download(
        self,
        file_url: str,
        filename: str,
        response_ok: t.Callable[[requests.Response], bool],
    ) -> str:
        assert self._cache_dir is not None
        os.makedirs(self._cache_dir, exist_ok=True)
        dest = os.path.join(self._cache_dir, filename)

        def check_response_for_download(r: requests.Response) -> bool:
            # if the response indicates a cache hit, treat it as valid
            # this ensures that we short-circuit any further evaluation immediately on
            # a hit
            if _cache_hit(dest, r):
                return True
            # we now know it's not a hit, so validate the content (forces download)
            return response_ok(r)

        response = _get_request(file_url, response_ok=check_response_for_download)
        # check to see if we have a file which matches the connection
        # only download if we do not (cache miss, vs hit)
        if not _cache_hit(dest, response):
            _atomic_write(dest, response.content)

        return dest

    @contextlib.contextmanager
    def open(
        self,
        file_url: str,
        filename: str,
        validate_response: t.Callable[[requests.Response], bool],
    ) -> t.Iterator[t.IO[bytes]]:
        if (not self._cache_dir) or self._disable_cache:
            yield io.BytesIO(
                _get_request(file_url, response_ok=validate_response).content
            )
        else:
            with open(
                self._download(file_url, filename, response_ok=validate_response), "rb"
            ) as fp:
                yield fp

    def bind(
        self,
        file_url: str,
        filename: str | None = None,
        validation_callback: t.Callable[[bytes], t.Any] | None = None,
    ) -> BoundCacheDownloader:
        return BoundCacheDownloader(
            file_url, self, filename=filename, validation_callback=validation_callback
        )


class BoundCacheDownloader:
    def __init__(
        self,
        file_url: str,
        downloader: CacheDownloader,
        *,
        filename: str | None = None,
        validation_callback: t.Callable[[bytes], t.Any] | None = None,
    ) -> None:
        self._file_url = file_url
        self._filename = filename or url_to_cache_filename(file_url)
        self._downloader = downloader
        self._validation_callback = validation_callback

    @contextlib.contextmanager
    def open(self) -> t.Iterator[t.IO[bytes]]:
        with self._downloader.open(
            self._file_url,
            self._filename,
            validate_response=self._validate_response,
        ) as fp:
            yield fp

    def _validate_response(self, response: requests.Response) -> bool:
        if not self._validation_callback:
            return True

        try:
            self._validation_callback(response.content)
            return True
        except ValueError:
            return False
