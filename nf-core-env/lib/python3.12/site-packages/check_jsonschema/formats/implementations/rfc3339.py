import re

# this regex is based on the one from the rfc3339-validator package
# credit to the original author
# original license:
#
#    MIT License
#
#    Copyright (c) 2019, Nicolas Aimetti
#
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.
#
# modifications have been made for additional corner cases and speed
RFC3339_REGEX = re.compile(
    r"""
    ^
    (?:\d{4})
    -
    (?:0[1-9]|1[0-2])
    -
    (?:[0-3]\d)
    (?:[Tt])
    (?:[01]\d|2[0123])
    :
    (?:[0-5]\d)
    :
    (?:[0-5]\d)
    # (optional) fractional seconds
    (?:[\.,]\d+)?
    # UTC or offset
    (?:
        [Zz]
        | [+-](?:[01]\d|2[0123]):[0-5]\d
    )
    $
""",
    re.VERBOSE | re.ASCII,
)


def validate(date_str: object) -> bool:
    """Validate a string as a RFC3339 date-time."""
    if not isinstance(date_str, str):
        return True
    if not RFC3339_REGEX.match(date_str):
        return False

    year, month, day = int(date_str[:4]), int(date_str[5:7]), int(date_str[8:10])

    if month in {4, 6, 9, 11}:
        max_day = 30
    elif month == 2:
        max_day = 29 if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else 28
    else:
        max_day = 31
    if not 1 <= day <= max_day:
        return False
    return True


if __name__ == "__main__":
    import timeit

    N = 100_000
    tests = (
        ("long_fracsec", "2018-12-31T23:59:59.8446519776713Z"),
        ("basic", "2018-12-31T23:59:59Z"),
        ("in_february", "2018-02-12T23:59:59Z"),
        ("in_february_invalid", "2018-02-29T23:59:59Z"),
        ("missing_t", "2018-12-31 23:59:59Z"),
        ("invalid_day", "2018-12-41T23:59:59Z"),
    )

    print("benchmarking")
    for name, val in tests:
        all_times = timeit.repeat(
            f"validate({val!r})", globals=globals(), repeat=3, number=N
        )
        print(f"{name} (valid={validate(val)}): {int(min(all_times) / N * 10**9)}ns")
