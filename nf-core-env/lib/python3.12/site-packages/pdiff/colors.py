import colorama

RESET = colorama.Style.RESET_ALL
FILE_HEADER = colorama.Fore.BLUE
HUNK_HEADER = colorama.Fore.CYAN

ADD_FG = colorama.Fore.GREEN
DELETE_FG = colorama.Fore.RED
CHANGE_FG = colorama.Fore.YELLOW

ADD_BG = colorama.Fore.BLACK + colorama.Back.GREEN
DELETE_BG = colorama.Fore.BLACK + colorama.Back.RED
CHANGE_BG = colorama.Fore.BLACK + colorama.Back.YELLOW

MARKERS_FG = {
    "\x00+": ADD_FG,
    "\x00-": DELETE_FG,
    "\x00^": CHANGE_FG,
    "\x01": RESET,
}

MARKERS_BG = {
    "\x00+": ADD_BG,
    "\x00-": DELETE_BG,
    "\x00^": CHANGE_BG,
    "\x01": RESET,
}


def colorize(text, color):
    return color + text + RESET
