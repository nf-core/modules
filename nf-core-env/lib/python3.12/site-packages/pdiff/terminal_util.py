import collections
import os

terminal_size = collections.namedtuple("terminal_size", ("columns", "lines"))


# `shutil.get_terminal_size` only checks file descriptor 1, which fails when
# e.g. stdout is piped to `less`, so instead we try file descriptors {0,1,2} via
# `os.get_terminal_size` (or in python 2, the raw `fcntl.ioctl` approach)
def get_terminal_size():
    terminal_size_fn = getattr(os, "get_terminal_size", _get_terminal_size)

    for fd in range(3):
        try:
            return terminal_size_fn(fd)
        except Exception:
            pass

    return terminal_size(columns=80, lines=24)


def _get_terminal_size(fd):
    import fcntl
    import struct
    import termios

    lines, columns = struct.unpack("hh", fcntl.ioctl(fd, termios.TIOCGWINSZ, "\x00\x00\x00\x00"))
    return terminal_size(columns=columns, lines=lines)
