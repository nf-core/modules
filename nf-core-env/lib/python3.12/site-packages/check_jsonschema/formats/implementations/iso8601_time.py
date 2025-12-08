import re

TIME_REGEX = re.compile(
    r"""
    ^
    (?:[01]\d|2[0123])
    :
    (?:[0-5]\d)
    :
    (?:[0-5]\d)
    # (optional) fractional seconds
    (?:(\.|,)\d+)?
    # UTC or offset
    (?:
        Z
        | z
        | [+-](?:[01]\d|2[0123]):[0-5]\d
    )
    $
""",
    re.VERBOSE | re.ASCII,
)


def validate(time_str: object) -> bool:
    if not isinstance(time_str, str):
        return False
    return bool(TIME_REGEX.match(time_str))


if __name__ == "__main__":
    import timeit

    N = 100_000
    tests = (
        ("basic", "23:59:59Z"),
        ("long_fracsec", "23:59:59.8446519776713Z"),
    )

    print("benchmarking")
    for name, val in tests:
        all_times = timeit.repeat(
            f"validate({val!r})", globals=globals(), repeat=3, number=N
        )
        print(f"{name} (valid={validate(val)}): {int(min(all_times) / N * 10**9)}ns")
