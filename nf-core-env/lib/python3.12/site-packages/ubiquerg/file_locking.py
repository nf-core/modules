import errno
import glob
import logging
import os

from contextlib import contextmanager
from .files import wait_for_lock, create_file_racefree
from signal import signal, SIGINT, SIGTERM
from typing import Union

PID = os.getpid()
READ = f"read-{PID}"
READ_GLOB = "read-*"
WRITE = "write"
UNIVERSAL = "universal"
LOCK_PREFIX = "lock"

_LOGGER = logging.getLogger(__name__)


class ThreeLocker(object):
    """
    A class to lock files for reading and writing.

    It uses a three-lock system, with separate read-lock, write-lock, and
    universal-lock (or lock-lock). The universal lock is used to lock the locks,
    to prevent race conditions between read and write locks. It allows multiple
    simultaneous readers, as long as there is no writer.
    It creates lock files in the same directory as the file to be locked.
    """

    def __init__(
        self, filepath: str, wait_max: int = 10, strict_ro_locks: bool = False
    ):
        self.wait_max = wait_max
        self.strict_ro_locks = strict_ro_locks
        self.set_file_path(filepath)
        self.locked = {READ: False, WRITE: False}

    @property
    def filepath(self) -> str:
        return self._filepath

    def set_file_path(self, filepath: str) -> str:
        if filepath:
            self._filepath = mkabs(filepath)
            self.lock_paths = make_all_lock_paths(self.filepath)
        else:
            self._filepath = None
            self.lock_paths = None

        return self._filepath

    def read_lock(self) -> bool:
        if not self.filepath:
            _LOGGER.warning("No filepath, no need to lock.")
            return True
        lock_path = self.lock_paths[READ]
        if not ensure_write_access(lock_path, self.strict_ro_locks):
            return True

        self.create_read_lock(self.filepath, self.wait_max)
        self.locked[READ] = True
        return True

    def write_lock(self) -> bool:
        if not self.filepath:
            _LOGGER.warning("No filepath, no need to lock.")
            return True

        lock_path = self.lock_paths[WRITE]
        if not ensure_write_access(lock_path, self.strict_ro_locks):
            # for writing, just fail anyway
            raise OSError(f"No write access to '{lock_path}'; can't lock file.")
        self.create_write_lock(self.filepath, self.wait_max)
        self.locked[READ] = True
        self.locked[WRITE] = True
        return True

    def read_unlock(self) -> bool:
        if not self.filepath:
            _LOGGER.warning("No filepath, no need to unlock.")
            return True
        _remove_lock(self.lock_paths[READ])
        self.locked[READ] = False
        return True

    def write_unlock(self) -> bool:
        if not self.filepath:
            _LOGGER.warning("No filepath, no need to unlock.")
            return True
        _remove_lock(self.lock_paths[WRITE])
        _remove_lock(self.lock_paths[READ])
        self.locked[WRITE] = False
        self.locked[READ] = False
        return True

    def create_read_lock(self, filepath: str = None, wait_max: int = None) -> None:
        """Securely create a read lock file.

        Args:
            filepath: path to a file to lock
            wait_max: max wait time if the file in question is already locked
        """
        filepath = filepath or self.filepath
        wait_max = wait_max or self.wait_max
        wait_for_lock(self.lock_paths[UNIVERSAL], wait_max)
        _create_lock(self.lock_paths[UNIVERSAL], filepath, wait_max)
        wait_for_lock(self.lock_paths[WRITE], wait_max)
        _create_lock(self.lock_paths[READ], filepath, wait_max)
        _remove_lock(self.lock_paths[UNIVERSAL])

    def create_write_lock(self, filepath: str = None, wait_max: int = None) -> None:
        """Securely create a write lock file.

        Args:
            filepath: path to a file to lock
            wait_max: max wait time if the file in question is already locked
        """
        filepath = filepath or self.filepath
        wait_max = wait_max or self.wait_max
        wait_for_lock(self.lock_paths[UNIVERSAL], wait_max)
        _create_lock(self.lock_paths[UNIVERSAL], filepath, wait_max)
        read_lock_paths = glob.glob(
            self.lock_paths[READ_GLOB]
        )  # must occur after universal lock is set
        all_lock_paths = read_lock_paths + [self.lock_paths[WRITE]]
        wait_for_locks(all_lock_paths, wait_max)
        _create_lock(self.lock_paths[READ], filepath, wait_max)
        _create_lock(self.lock_paths[WRITE], filepath, wait_max)
        _remove_lock(self.lock_paths[UNIVERSAL])

    def _interrupt_handler(self, signal_received, frame):
        if signal_received == SIGINT:
            _LOGGER.warning(f"Received SIGINT, unlocking file and exiting...")
            self.write_unlock()
            self.read_unlock()
            raise SystemExit
        if signal_received == SIGTERM:
            _LOGGER.warning(f"Received SIGTERM, unlocking file and exiting...")
            self.__exit__(None, None, None)
            self.write_unlock()
            self.read_unlock()
            raise SystemExit

    def __repr__(self) -> str:
        settings_dict = {
            "filepath": self.filepath,
            "wait_max": self.wait_max,
            "locked": self.locked,
            "strict_ro_locks": self.strict_ro_locks,
        }

        return f"{type(self).__name__}({settings_dict})"

    def __del__(self) -> None:
        if self.filepath:
            if self.locked[WRITE]:
                self.write_unlock()
            if self.locked[READ]:
                self.read_unlock()


def ensure_locked(type: str = WRITE):  # decorator factory
    """Decorator to apply to functions to make sure they only happen when locked."""

    def decorator(func):
        def inner_func(self, *args, **kwargs):
            if not self.locker:
                raise OSError("File not lockable. File locker not provided.")
            if not self.locker.locked[type]:
                raise OSError(
                    f"This function must use a context manager to {type}-lock the file"
                )

            return func(self, *args, **kwargs)

        return inner_func

    return decorator


@contextmanager
def read_lock(obj: Union[str, object]) -> object:
    """Read-lock a filepath or object with locker attribute.

    Args:
        obj: filepath string or object with locker attribute

    Yields:
        object: the locked object
    """
    if type(obj) == str:
        locker = ThreeLocker(obj)
    elif hasattr(obj, "locker"):
        locker = obj.locker
    else:
        raise AttributeError(f"Cannot lock: {obj}.")

    # handle a premature Ctrl+C exit from this context manager
    try:
        signal(SIGTERM, locker._interrupt_handler)
        signal(SIGINT, locker._interrupt_handler)
        # If this is run in a thread, the signal module is not available and raises an exception.
        # ValueError: signal only works in main thread of the main interpreter
        # That's fine; in this case, we don't need to handle signals anyway.
    except ValueError as e:
        _LOGGER.error(f"Failed to set interrupt handler: {e}")

    locker.read_lock()

    try:
        yield obj
    finally:
        locker.read_unlock()


@contextmanager
def write_lock(obj: Union[str, object]) -> object:
    """Write-lock file path or object with locker attribute.

    Args:
        obj: filepath string or object with locker attribute

    Yields:
        object: the locked object
    """
    if type(obj) == str:
        locker = ThreeLocker(obj)
    elif hasattr(obj, "locker"):
        locker = obj.locker
    else:
        raise AttributeError(f"Cannot lock: {obj}.")

    # handle a premature Ctrl+C exit from this context manager
    try:
        signal(SIGTERM, locker._interrupt_handler)
        signal(SIGINT, locker._interrupt_handler)
    except ValueError as e:
        _LOGGER.error(f"Failed to set interrupt handler: {e}")

    locker.write_lock()
    try:
        yield obj
    finally:
        locker.write_unlock()


def locked_read_file(filepath, create_file: bool = False) -> str:
    """Read a file contents into memory after locking the file.

    This will prevent other ThreeLocker-protected processes from writing to the
    file while it is being read.

    Args:
        filepath: path to the file that should be read
        create_file: whether to create the file if it doesn't exist

    Returns:
        str: file contents
    """
    if os.path.exists(filepath):
        with read_lock(filepath), open(filepath, "r") as file:
            file_contents = file.read()
    elif create_file:
        _LOGGER.info("File does not exist, but create_file is true. Creating...")
        file_contents = ""
        create_file_racefree(filepath)
    else:
        raise FileNotFoundError(f"No such file: {filepath}")
    return file_contents


def wait_for_locks(lock_paths: Union[list, str], wait_max: int = 10):
    """Wait for lock files to be removed.

    Args:
        lock_paths: path to a file to lock
        wait_max: max wait time if the file in question is already locked
    """
    if not isinstance(lock_paths, list):
        lock_paths = [lock_paths]
    for lock_path in lock_paths:
        wait_for_lock(lock_path, wait_max)


def ensure_write_access(lock_path, strict_ro_locks=False):
    if os.access(os.path.dirname(lock_path), os.W_OK):  # write access
        return True
    else:
        if strict_ro_locks:  # fail; no write access, strict mode
            raise OSError(f"No write access to '{lock_path}'; can't lock file.")
        else:  # warn; no write access, non-strict mode
            _LOGGER.warning(f"No write access to '{lock_path}'; can't lock file.")
            return True


def _create_lock(lock_path, filepath, wait_max) -> None:
    # _LOGGER.debug(f"Creating lock for {filepath} at {lock_path}")
    _LOGGER.debug(f"Creating lock at {os.path.basename(lock_path)}")
    try:
        create_file_racefree(lock_path)
    except FileNotFoundError:
        parent_dir = os.path.dirname(filepath)
        os.makedirs(parent_dir)
        _create_lock(lock_path, filepath, wait_max)
    except Exception as e:
        if e.errno == errno.EEXIST:
            # Rare case: file already exists;
            # the lock has been created in the split second since the
            # last lock existence check,
            # wait for the lock file to be gone, but no longer than
            # `wait_max`.
            _LOGGER.info(
                "The lock has been created in the split second since the "
                "last lock existence check. Waiting"
            )
            wait_for_lock(lock_path, wait_max)
            _create_lock(lock_path, filepath, wait_max)
        else:
            raise e


def _remove_lock(lock_path) -> bool:
    """Remove lock.

    Args:
        lock_path: path to the lock file to remove. Not the path to the actual file!

    Returns:
        bool: whether the lock was found and removed
    """
    _LOGGER.debug(f"Removing lock at {os.path.basename(lock_path)}")
    if os.path.exists(lock_path):
        os.remove(lock_path)
        return True
    return False


def make_all_lock_paths(filepath):
    """
    Create a collection of paths to lock files with given name as base.
    """
    lock_paths = {}
    for type in [READ, WRITE, UNIVERSAL, READ_GLOB]:
        prefix = f"{LOCK_PREFIX}-{type}-" if type else LOCK_PREFIX
        base, name = os.path.split(filepath)
        lock_name = name if name.startswith(prefix) else prefix + name
        lock_paths[type] = lock_name if not base else os.path.join(base, lock_name)
    return lock_paths


def mkabs(path, reldir=None):
    """Make sure a path is absolute.

    If not already absolute, it's made absolute relative to a given directory.
    Also expands ~ and environment variables for kicks.

    Args:
        path: Path to make absolute
        reldir: Relative directory to make path absolute from if it's not already absolute

    Returns:
        str: Absolute path
    """

    def xpand(path):
        return os.path.expandvars(os.path.expanduser(path))

    if os.path.isabs(xpand(path)):
        return xpand(path)

    if not reldir:
        return os.path.abspath(xpand(path))

    return os.path.join(xpand(reldir), xpand(path))


# class OneLocker(object):
# def lock(self, type=READ):
#     if not self.filepath:
#         _LOGGER.warning("No filepath, no need to lock.")
#         return True

#     # Check for permissions to write a lock file
#     lock_path = self.lock_paths[type]

#     if not os.access(os.path.dirname(lock_path), os.W_OK):
#         if self.strict_ro_locks:
#             raise OSError(f"No write access to '{lock_path}'; can't lock file.")
#         else:
#             _LOGGER.warning(f"No write access to '{lock_path}'; can't lock file.")
#             self.locked[type] = True
#             return True

#     create_lock(self.filepath, self.wait_max)
#     self.locked[type] = True
#     return True

# def unlock(self, type=READ):
#     if not self.filepath:
#         _LOGGER.warning("No filepath, no need to unlock.")
#         return True

#     # Check for permissions to write a lock file
#     lock_path = make_lock_path(self.filepath)
#     if not os.access(os.path.dirname(lock_path), os.W_OK):
#         if self.strict_ro_locks:
#             raise OSError(f"No write access to '{lock_path}' can't lock file.")
#         else:
#             _LOGGER.warning(f"No write access to '{lock_path}' can't lock file.")
#             self.locked = False
#             return True

#     remove_lock(self.filepath)
#     self.locked = False
#     return True

# def __del__(self):
#     if self.filepath and self.locked:
#         self.unlock()
