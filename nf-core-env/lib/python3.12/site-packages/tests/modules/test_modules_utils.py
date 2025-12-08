import nf_core.modules.modules_utils

from ..test_modules import TestModules


class TestModulesUtils(TestModules):
    def test_get_installed_modules(self):
        """Test getting installed modules"""
        _, nfcore_modules = nf_core.modules.modules_utils.get_installed_modules(self.nfcore_modules)
        assert len(nfcore_modules) == 1
        assert nfcore_modules[0].component_name == "bpipe/test"

    def test_get_installed_modules_with_files(self):
        """Test getting installed modules. When a module contains a file in its directory, it shouldn't be picked up as a tool/subtool"""
        # Create a file in the module directory
        with open(self.nfcore_modules / "modules" / "nf-core" / "bpipe" / "test_file.txt", "w") as fh:
            fh.write("test")

        _, nfcore_modules = nf_core.modules.modules_utils.get_installed_modules(self.nfcore_modules)
        assert len(nfcore_modules) == 1

    def test_filter_modules_by_name_exact_match(self):
        """Test filtering modules by name with an exact match"""
        # install bpipe/test
        _, nfcore_modules = nf_core.modules.modules_utils.get_installed_modules(self.nfcore_modules)

        # Test exact match
        filtered = nf_core.modules.modules_utils.filter_modules_by_name(nfcore_modules, "bpipe/test")
        assert len(filtered) == 1
        assert filtered[0].component_name == "bpipe/test"

    def test_filter_modules_by_name_tool_family(self):
        """Test filtering modules by name to get all subtools of a super-tool"""
        # Create some mock samtools subtools in the modules directory
        samtools_dir = self.nfcore_modules / "modules" / "nf-core" / "samtools"

        for subtool in ["view", "sort", "index"]:
            subtool_dir = samtools_dir / subtool
            subtool_dir.mkdir(parents=True, exist_ok=True)
            (subtool_dir / "main.nf").touch()

        # Get the modules
        _, nfcore_modules = nf_core.modules.modules_utils.get_installed_modules(self.nfcore_modules)

        # Test filtering by tool family (super-tool)
        filtered = nf_core.modules.modules_utils.filter_modules_by_name(nfcore_modules, "samtools")

        assert set(m.component_name for m in filtered) == {"samtools/view", "samtools/sort", "samtools/index"}

    def test_filter_modules_by_name_exact_match_preferred(self):
        """Test that exact matches are preferred over prefix matches"""
        # Create a samtools super-tool and its subtools
        samtools_dir = self.nfcore_modules / "modules" / "nf-core" / "samtools"
        samtools_dir.mkdir(parents=True, exist_ok=True)
        (samtools_dir / "main.nf").touch()

        # Create subtools
        for subtool in ["view", "sort"]:
            subtool_dir = samtools_dir / subtool
            subtool_dir.mkdir(parents=True, exist_ok=True)
            (subtool_dir / "main.nf").touch()

        # Get the modules
        _, nfcore_modules = nf_core.modules.modules_utils.get_installed_modules(self.nfcore_modules)

        # Test that exact match is returned when it exists
        filtered = nf_core.modules.modules_utils.filter_modules_by_name(nfcore_modules, "samtools")
        assert len(filtered) == 1
        assert filtered[0].component_name == "samtools"

    def test_filter_modules_by_name_no_match(self):
        """Test filtering modules by name with no matches"""
        _, nfcore_modules = nf_core.modules.modules_utils.get_installed_modules(self.nfcore_modules)

        # Test no match
        filtered = nf_core.modules.modules_utils.filter_modules_by_name(nfcore_modules, "nonexistent")
        assert len(filtered) == 0

    def test_filter_modules_by_name_empty_list(self):
        """Test filtering an empty list of modules"""
        modules = []

        filtered = nf_core.modules.modules_utils.filter_modules_by_name(modules, "fastqc")
        assert len(filtered) == 0
