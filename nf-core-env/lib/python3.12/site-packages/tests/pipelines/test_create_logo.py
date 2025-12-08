"""Test covering the create-logo command."""

from pathlib import Path

import nf_core.pipelines.create_logo

from ..test_pipelines import TestPipelines


class TestCreateLogo(TestPipelines):
    def test_create_logo_png(self):
        """Test that the create-logo command works for PNGs"""

        # Create a logo
        logo_fn = nf_core.pipelines.create_logo.create_logo("pipes", self.pipeline_dir)
        # Check that the file exists
        self.assertTrue(logo_fn.is_file())
        # Check that the file is a PNG
        self.assertTrue(logo_fn.suffix == ".png")
        # Check that the file is the right size
        fixture_fn = Path(__file__).parent.parent / "fixtures" / "create_logo.png"
        # allow some flexibility in the file size
        self.assertTrue(int(logo_fn.stat().st_size / 1000) == int(fixture_fn.stat().st_size / 1000))

    def test_create_logo_png_dark(self):
        """Test that the create-logo command works for dark PNGs"""

        # Create a logo
        logo_fn = nf_core.pipelines.create_logo.create_logo("pipes", self.pipeline_dir, theme="dark")
        # Check that the file exists
        self.assertTrue(logo_fn.is_file())
        # Check that the file is a PNG
        self.assertTrue(logo_fn.suffix == ".png")
        # Check that the file is the right size
        fixture_fn = Path(__file__).parent.parent / "fixtures" / "create_logo_dark.png"
        # allow some flexibility in the file size
        self.assertTrue(int(logo_fn.stat().st_size / 1000) == int(fixture_fn.stat().st_size / 1000))

    def test_create_log_png_width(self):
        """Test that the create-logo command works for PNGs with a custom width"""

        # Create a logo
        logo_fn = nf_core.pipelines.create_logo.create_logo("pipes", self.pipeline_dir, width=100)
        # Check that the file exists
        self.assertTrue(logo_fn.is_file())
        # Check that the file is a PNG
        self.assertTrue(logo_fn.suffix == ".png")
        # Check that the file is the right size
        fixture_fn = Path(__file__).parent.parent / "fixtures" / "create_logo_width100.png"
        # allow some flexibility in the file size
        self.assertTrue(int(logo_fn.stat().st_size / 100) == int(fixture_fn.stat().st_size / 100))

    def test_create_logo_twice(self):
        """Test that the create-logo command returns an info message when run twice"""

        # Create a logo
        logo_fn = nf_core.pipelines.create_logo.create_logo("pipes", self.pipeline_dir)
        # Check that the file exists
        self.assertTrue(logo_fn.is_file())
        # Create the logo again and capture the log output
        with self.assertLogs(level="INFO") as log:
            nf_core.pipelines.create_logo.create_logo("pipes", self.pipeline_dir)
            # Check that the log message is correct
            self.assertIn("Logo already exists", log.output[0])

    def test_create_logo_without_text_fail(self):
        """Test that the create-logo command fails without text"""

        # Create a logo
        with self.assertRaises(UserWarning):
            nf_core.pipelines.create_logo.create_logo("", self.pipeline_dir)

    def test_create_logo_with_filename(self):
        """Test that the create-logo command works with a custom filename"""

        # Create a logo
        logo_fn = nf_core.pipelines.create_logo.create_logo(
            "pipes", Path(self.pipeline_dir / "custom_dir"), filename="custom"
        )
        # Check that the file exists
        self.assertTrue(logo_fn.is_file())
        # Check that the parent directory name
        self.assertTrue(logo_fn.parent.name == "custom_dir")
        # Check that the file has correct name
        self.assertTrue(logo_fn.name == "custom.png")

    def test_create_logo_svg(self):
        """Test that the create-logo command works for SVGs"""

        # Create a logo
        logo_fn = nf_core.pipelines.create_logo.create_logo("pipes", self.pipeline_dir, format="svg")
        # Check that the file exists
        self.assertTrue(logo_fn.is_file())
        # Check that the file is a SVG
        self.assertTrue(logo_fn.suffix == ".svg")
        # Check that the svg contains the correct text
        with open(logo_fn) as fh:
            svg = fh.read()
        self.assertIn("pipes", svg)
        # check that it is the light theme
        self.assertIn("#050505", svg)

    def test_create_logo_svg_dark(self):
        """Test that the create-logo command works for svgs and dark theme"""

        # Create a logo
        logo_fn = nf_core.pipelines.create_logo.create_logo("pipes", self.pipeline_dir, format="svg", theme="dark")
        # Check that the file exists
        self.assertTrue(logo_fn.is_file())
        # Check that the file is a SVG
        self.assertTrue(logo_fn.suffix == ".svg")
        # Check that the svg contains the correct text
        with open(logo_fn) as fh:
            svg = fh.read()
        self.assertIn("pipes", svg)
        # check that it is the dark theme
        self.assertIn("#fafafa", svg)
