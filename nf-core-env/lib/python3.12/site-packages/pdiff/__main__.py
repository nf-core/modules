import os
import sys

from . import argument_parser
from . import diff_formatter


def main():
    parser = argument_parser.get_parser()
    args = parser.parse_args()

    for filename in (args.left_filename, args.right_filename):
        if not os.path.exists(filename):
            sys.stderr.write("error: file does not exist: %s\n" % filename)
            return 1
        elif os.path.isdir(filename):
            sys.stderr.write("error: path is a directory: %s\n" % filename)
            return 1

    formatter = diff_formatter.DiffFormatter(
        args.left_filename,
        args.right_filename,
        args.context,
        args.width,
        args.tab_size,
        args.signs,
        args.line_numbers,
        args.background,
    )

    for line in formatter.get_lines():
        sys.stdout.write(line)

    return 0


if __name__ == "__main__":
    sys.exit(main())
