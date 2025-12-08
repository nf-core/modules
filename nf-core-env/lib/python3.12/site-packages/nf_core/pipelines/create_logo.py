import logging
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

import nf_core
from nf_core.utils import NFCORE_CACHE_DIR

log = logging.getLogger(__name__)


def create_logo(
    text: str,
    directory: Path | str,
    filename: str = "",
    theme: str = "light",
    width: int = 2300,
    format: str = "png",
    force: bool = False,
) -> Path:
    """Create a logo for a pipeline."""
    if not text:
        raise UserWarning("Please provide the name of the text to put on the logo.")
    directory = Path(directory)
    if not directory.is_dir():
        log.debug(f"Creating directory {directory}")
        directory.mkdir(parents=True, exist_ok=True)
    assets = Path(nf_core.__file__).parent / "assets/logo"

    if format == "svg":
        template_fn = "placeholder_logo.svg"

        if width != 2300:
            log.warning("SVG format does not support resizing. Setting width to 2300px.")

        # replace the placeholder text with the pipeline name
        with open(assets / template_fn) as fh:
            svg = fh.read().replace("PLACEHOLDER", text)
        if theme == "dark":
            svg = svg.replace("#050505", "#fafafa")

        # save the svg
        logo_filename = f"nf-core-{text}_logo_{theme}.svg" if not filename else filename
        logo_filename = f"{logo_filename}.svg" if not logo_filename.lower().endswith(".svg") else logo_filename
        logo_path = Path(directory, logo_filename)
        with open(logo_path, "w") as fh:
            fh.write(svg)

    else:
        logo_filename = f"nf-core-{text}_logo_{theme}.png" if not filename else filename
        logo_filename = f"{logo_filename}.png" if not logo_filename.lower().endswith(".png") else logo_filename
        cache_name = f"nf-core-{text}_logo_{theme}_{width}.png"
        logo_path = Path(directory, logo_filename)

        # Check if we haven't already created this logo
        if logo_path.is_file() and not force:
            log.info(f"Logo already exists at: {logo_path}. Use `--force` to overwrite.")
            return logo_path
        # cache file
        cache_path = Path(NFCORE_CACHE_DIR, "logo", cache_name)
        img: Image.Image | None = None
        if cache_path.is_file():
            log.debug(f"Logo already exists in cache at: {cache_path}. Reusing this file.")
            img = Image.open(cache_path)
        if img is None:
            log.debug(f"Creating logo for {text}")

            # make sure the figure fits the text
            font_path = assets / "MavenPro-Bold.ttf"
            log.debug(f"Using font: {str(font_path)}")
            font = ImageFont.truetype(str(font_path), 400)
            text_length = font.getmask(text).getbbox()[2]  # get the width of the text based on the font

            max_width = max(
                2300, text_length + len(text) * 20
            )  # need to add some more space to the text length to make sure it fits

            template_fn = "nf-core-repo-logo-base-lightbg.png"
            if theme == "dark":
                template_fn = "nf-core-repo-logo-base-darkbg.png"

            template_path = assets / template_fn
            img = Image.open(template_path)
            if img is None:
                raise RuntimeError("Failed to create logo image")
            # get the height of the template image
            height = img.size[1]

            # Draw text
            draw = ImageDraw.Draw(img)
            color = theme == "dark" and (250, 250, 250) or (5, 5, 5)
            draw.text((110, 465), text, color, font=font)

            if img is not None:
                # Crop to max width
                img = img.crop((0, 0, max_width, height))

                # Resize
                img = img.resize((width, int((width / max_width) * height)))
            else:
                log.error("Failed to create logo, no image object created.")

            # Save to cache
            Path(cache_path.parent).mkdir(parents=True, exist_ok=True)
            log.debug(f"Saving logo to cache: {cache_path}")
            if img is not None:
                img.save(cache_path, "PNG")
        # Save
        if img is not None:
            img.save(logo_path, "PNG")
        else:
            raise RuntimeError("Failed to create logo image")

    log.debug(f"Saved logo to: '{logo_path}'")

    # Return the logo
    return logo_path
