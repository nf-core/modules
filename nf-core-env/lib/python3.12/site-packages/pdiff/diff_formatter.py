import difflib
import itertools
import re

from . import colors


class DiffFormatter(object):
    LEFT_HEADER_PREFIX = "--- "
    RIGHT_HEADER_PREFIX = "+++ "
    LEFT_HUNK_HEADER_TEMPLATE = "@@ -%d,%d @@"
    RIGHT_HUNK_HEADER_TEMPLATE = "@@ +%d,%d @@"

    def __init__(
        self,
        left_filename,
        right_filename,
        context,
        width,
        tab_size,
        signs,
        line_numbers,
        background,
    ):
        self.left_filename = left_filename
        self.right_filename = right_filename
        self.context = context
        self.width = width
        self.tab_size = tab_size
        self.signs = signs
        self.line_numbers = line_numbers
        self.background = background

        self.half_width = self.width // 2 - 1
        self.empty_half = " " * self.half_width
        self.tab_spaces = " " * self.tab_size
        self.marker_to_colors = colors.MARKERS_BG if self.background else colors.MARKERS_FG

    def get_lines(self):
        # `readlines` because `difflib._mdiff` can't operate on a generator
        with open(self.left_filename, "r") as left_file:
            left_lines = left_file.readlines()
        with open(self.right_filename, "r") as right_file:
            right_lines = right_file.readlines()

        yield self._format_file_header()

        # Determine line number column widths; note: in some cases this is
        # unnecessarily wide (e.g. if changes only in beginning of a long file)
        left_ln_width = len(str(len(left_lines)))
        right_ln_width = len(str(len(right_lines)))

        mdiff = difflib._mdiff(left_lines, right_lines, context=self.context)

        # `mdiff` context separators don't have the metadata necessary to generate
        # git-diff-like hunk headers (`@@ -%d,%d @@` and `@@ +%d,%d @@`), so we
        # partition `mdiff` into hunks and process each one separately
        for hunk in self._get_hunks(mdiff):
            for line in self._format_hunk(hunk, left_ln_width, right_ln_width):
                yield line

    def _format_file_header(self):
        left_header = colors.colorize(
            DiffFormatter.LEFT_HEADER_PREFIX + self.left_filename, colors.FILE_HEADER
        )
        right_header = colors.colorize(
            DiffFormatter.RIGHT_HEADER_PREFIX + self.right_filename, colors.FILE_HEADER
        )
        return self._format_line(left_header, right_header)

    def _get_hunks(self, mdiff):
        hunk = []

        for mdiff_tuple in mdiff:
            if mdiff_tuple[2] is None:
                if hunk:
                    yield hunk
                    hunk = []
            else:
                hunk.append(mdiff_tuple)

        # Don't forget the last hunk, which isn't followed by a context separator
        if hunk:
            yield hunk

    def _format_hunk(self, hunk, left_ln_width, right_ln_width):
        yield self._format_hunk_header(hunk)

        for (left_num, left_half), (right_num, right_half), has_changes in hunk:
            left_half = left_half.replace("\n", "").replace("\t", self.tab_spaces)
            right_half = right_half.replace("\n", "").replace("\t", self.tab_spaces)

            if has_changes:
                for marker, color in self.marker_to_colors.items():
                    left_half = left_half.replace(marker, color)
                    right_half = right_half.replace(marker, color)

                # Use background color for whitespace-only diffs even if no --background
                if left_half and right_half and not self.background:
                    left_half = self._highlight_whitespace_background(left_half)
                    right_half = self._highlight_whitespace_background(right_half)

            left_sign, right_sign = self._format_signs(left_half, right_half, has_changes)

            if self.line_numbers:
                left_half = str(left_num).rjust(left_ln_width) + " " + left_half
                right_half = str(right_num).rjust(right_ln_width) + " " + right_half

            yield self._format_line(left_sign + left_half, right_sign + right_half)

    def _format_hunk_header(self, hunk):
        left_nums = [left_num for (left_num, _), _, _ in hunk if left_num != ""]
        right_nums = [right_num for _, (right_num, _), _ in hunk if right_num != ""]
        left_start = left_nums[0] if left_nums else 0
        right_start = right_nums[0] if right_nums else 0
        left_header = colors.colorize(
            DiffFormatter.LEFT_HUNK_HEADER_TEMPLATE % (left_start, len(left_nums)),
            colors.HUNK_HEADER,
        )
        right_header = colors.colorize(
            DiffFormatter.RIGHT_HUNK_HEADER_TEMPLATE % (right_start, len(right_nums)),
            colors.HUNK_HEADER,
        )
        return self._format_line(left_header, right_header)

    def _highlight_whitespace_background(self, line):
        # For whitespace-only diffs, change `\x1b[30m` to `\x1b[40m` (background)
        return re.sub("(?<=\x1b\\[)3(?=[0-7]m\\s+\x1b\\[0m)", "4", line)

    def _format_signs(self, left_half, right_half, has_changes):
        if not self.signs:
            return "", ""
        if not has_changes:
            return "  ", "  "
        elif not left_half:
            return "  ", colors.colorize("+ ", colors.ADD_FG)
        elif not right_half:
            return colors.colorize("- ", colors.DELETE_FG), "  "
        else:
            return (colors.colorize("! ", colors.CHANGE_FG),) * 2

    def _format_line(self, left_half, right_half):
        left_half_lines = self._format_half_lines(left_half)
        right_half_lines = self._format_half_lines(right_half)

        return (
            "\n".join(
                left_half + " " + right_half
                for left_half, right_half in itertools.zip_longest(
                    left_half_lines, right_half_lines, fillvalue=self.empty_half
                )
            )
            + "\n"
        )

    def _format_half_lines(self, half_line):
        # Split `'ab\x1b[31mcd\x1b[0m'` into `['ab', '\x1b[31m', 'cd', '\x1b[0m']`
        parts = re.split("(\x1b\\[(?:0|[34][0-7])m)", half_line)
        visible_len = sum(len(part) for part in parts if not part.startswith("\x1b"))

        if visible_len <= self.half_width:
            half_lines = [half_line]
        else:
            half_lines = [""]
            visible_len = 0
            last_color = colors.RESET

            # We can't simply wrap `half_line` at every `half_width` chars due to
            # color codes, so we manually assemble `half_lines` from `parts`
            for part in parts:
                if part.startswith("\x1b"):
                    half_lines[-1] += part
                    last_color += part
                elif visible_len + len(part) <= self.half_width:
                    half_lines[-1] += part
                    visible_len += len(part)
                else:
                    first_offset = self.half_width - visible_len
                    half_lines[-1] += part[:first_offset] + colors.RESET

                    for offset in range(first_offset, len(part), self.half_width):
                        wrapped_half_line = part[offset : offset + self.half_width]
                        half_lines.append(last_color + wrapped_half_line + colors.RESET)
                        visible_len = len(wrapped_half_line)

        # Always right pad the last half line
        pad_len = self.half_width - visible_len
        half_lines[-1] += " " * pad_len
        return half_lines
