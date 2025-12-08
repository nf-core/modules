"""Generate HTML reports"""

import csv
import glob
import os
import shutil
import sys
from copy import deepcopy
from datetime import timedelta
from json import dumps
from logging import getLogger
from typing import List

import jinja2
import pandas as _pd
import yaml
from eido import read_schema
from peppy.const import AMENDMENTS_KEY
from ubiquerg import mkabs

from ._version import __version__
from .const import (
    BUTTON_APPEARANCE_BY_FLAG,
    FILE_KEY,
    NO_DATA_PLACEHOLDER,
    OBJECT_TYPES,
    OUTPUT_DIR,
    PIPELINE_NAME,
    PIPELINE_TYPE,
    PKG_NAME,
    PROFILE_COLNAMES,
    PROJECT_NAME,
    STATUS_FILE_DIR,
    TEMPLATES_DIRNAME,
)
from .helpers import make_subdirectories

_LOGGER = getLogger(PKG_NAME)


class HTMLReportBuilder(object):
    """Generate HTML summary report for project/samples"""

    def __init__(self, prj, portable=False):
        """
        The Project defines the instance.

        :param PipestatManager prj: Project with which to work/operate on
        :param bool portable: moves figures and report files to directory for easy sharing
        """

        self.prj = prj
        self.jinja_env = get_jinja_env()

        self.portable = portable

        results_file_path = getattr(self.prj.backend, "results_file_path", None)
        config_path = self.prj.cfg.get("config_path", None)
        output_dir = self.prj.cfg.get(OUTPUT_DIR, None)
        self.output_dir = output_dir or results_file_path or config_path

        if os.path.isdir(self.output_dir):
            pass
        else:
            self.output_dir = os.path.dirname(self.output_dir)

        if not self.portable:
            self.reports_dir = os.path.join(self.output_dir, "reports")
        else:
            _LOGGER.info(f"Portable flag set: {self.portable}")
            self.reports_dir = os.path.join(self.output_dir, "portable_reports")

        _LOGGER.debug(f"Reports dir: {self.reports_dir}")

    def __call__(
        self, pipeline_name, project_index_html=None, amendment=None, looper_samples=None
    ):
        """
        Generate HTML report.

        :param str pipeline_name: ID of the pipeline to generate the report for
        :param Iterable[str] amendment: name indicating amendment to use, optional
        :return str: path to the index page of the generated HTML report
        """
        # Generate HTML report
        self.pipeline_name = pipeline_name
        self.amendment = amendment
        self.amendments_str = "_".join(self.amendment) if self.amendment else ""
        self.pipeline_reports = os.path.join(
            self.reports_dir,
            (
                f"{self.pipeline_name}_{self.amendments_str}"
                if self.amendments_str
                else self.pipeline_name
            ),
        )
        self.prj_index_html_path = project_index_html
        self.index_html_path = os.path.join(self.pipeline_reports, "index.html")
        self.looper_samples = looper_samples or None

        self.schema = deepcopy(self.prj.result_schemas)
        for k, v in self.schema.items():
            if "object_type" in self.schema[k]:
                self.schema[k]["type"] = self.schema[k]["object_type"]

        navbar = self.create_navbar(
            navbar_links=self.create_navbar_links(
                wd=self.pipeline_reports,
                project_index_html_relpath=(
                    os.path.relpath(self.prj_index_html_path, self.pipeline_reports)
                    if self.prj_index_html_path
                    else None
                ),
            ),
            index_html_relpath=os.path.relpath(self.index_html_path, self.pipeline_reports),
        )
        self.create_index_html(navbar, self.create_footer())
        return self.index_html_path

    def _reset_pipeline_type(self):
        """
        The report logic will set the pipeline type when multi results is used. It must be reset or it causes issues.
        """
        self.prj.backend.pipeline_type = self.prj.pipeline_type

    def _create_copy_for_porting(self, parent_path: str, record_identifier: str) -> str:
        """
        Helper function that assists with copying images (pdfs)
        Original images stay in their original location.
        This does NOT move the original data.

        :param str parentpath: path of existing file
        :param str newpath: desired destination
        :return str newpath if successful else returns parentpath
        """
        destination_dir = os.path.join(self.pipeline_reports, "resources")
        full_name = os.path.basename(parent_path)

        destination_file = os.path.join(
            os.path.join(destination_dir, record_identifier), full_name
        )

        os.makedirs(os.path.dirname(destination_file), exist_ok=True)

        try:
            shutil.copyfile(parent_path, destination_file)
        except Exception as e:
            # if any exception occurs, simply return original path
            _LOGGER.warning(e)
            return parent_path

        return destination_file

    def create_object_parent_html(self, navbar, footer):
        """
        Generates a page listing all the project objects with links
        to individual object pages

        :param str navbar: HTML to be included as the navbar in the main summary page
        :param str footer: HTML to be included as the footer
        :return str: Rendered parent objects HTML file
        """
        if not os.path.exists(self.pipeline_reports):
            os.makedirs(self.pipeline_reports)
        pages = []
        labels = []
        obj_result_ids = self.get_nonhighlighted_results(OBJECT_TYPES)

        for key in obj_result_ids:
            desc = (
                self.prj.result_schemas[key]["description"]
                if "description" in self.prj.result_schemas[key]
                else ""
            )
            labels.append(f"<b>{key.replace('_', ' ')}</b>: {desc}")
            page_path = os.path.join(
                self.pipeline_reports, f"{key}.html".replace(" ", "_").lower()
            )
            pages.append(os.path.relpath(page_path, self.pipeline_reports))

        template_vars = dict(
            navbar=navbar,
            footer=footer,
            labels=labels,
            pages=pages,
            header="Objects",
            pipeline_name=self.pipeline_name,
        )
        _LOGGER.debug(f"object navbar_list_parent.html | template_vars:" f"\n{template_vars}")
        return render_jinja_template("navbar_list_parent.html", self.jinja_env, template_vars)

    def create_sample_parent_html(self, navbar, footer):
        """
        Generates a page listing all the project samples with links
        to individual sample pages
        :param str navbar: HTML to be included as the navbar in the main summary page
        :param str footer: HTML to be included as the footer
        :return str: Rendered parent samples HTML file
        """
        if not os.path.exists(self.pipeline_reports):
            os.makedirs(self.pipeline_reports)
        pages = []
        labels = []

        if self.prj.cfg["multi_result_files"] is True:
            pipeline_types = ["sample", "project"]
        else:
            pipeline_types = [self.prj.backend.pipeline_type]

        for pipeline_type in pipeline_types:
            self.prj.backend.pipeline_type = pipeline_type
            for sample in self.prj.backend.select_records()["records"]:
                sample_name = sample["record_identifier"]
                sample_dir = self.pipeline_reports

                # Confirm sample directory exists, then build page
                if os.path.exists(sample_dir):
                    page_path = os.path.join(
                        self.pipeline_reports,
                        f"{sample_name}.html".replace(" ", "_").lower(),
                    )
                    page_relpath = os.path.relpath(page_path, self.pipeline_reports)
                    pages.append(page_relpath)
                    labels.append(sample_name)

        self._reset_pipeline_type()
        template_vars = dict(
            navbar=navbar,
            footer=footer,
            labels=labels,
            pages=pages,
            header="Records",
            pipeline_name=self.pipeline_name,
        )
        _LOGGER.debug(f"sample navbar_list_parent.html | template_vars:" f"\n{template_vars}")
        return render_jinja_template("navbar_list_parent.html", self.jinja_env, template_vars)

    def create_navbar(self, navbar_links, index_html_relpath):
        """
        Creates the navbar using the provided links

        :param str navbar_links: HTML list of links to be inserted into a navbar
        :return str: navbar HTML
        """
        template_vars = dict(navbar_links=navbar_links, index_html=index_html_relpath)
        return render_jinja_template("navbar.html", self.jinja_env, template_vars)

    def create_footer(self):
        """
        Renders the footer from the templates directory

        :return str: footer HTML
        """
        return render_jinja_template("footer.html", self.jinja_env, dict(version=__version__))

    def create_navbar_links(self, wd=None, context=None, project_index_html_relpath=None):
        """
        Return a string containing the navbar prebuilt html.

        Generates links to each page relative to the directory of interest
        (wd arg) or uses the provided context to create the paths (context arg)

        :param path wd: the working directory of the current HTML page being
            generated, enables navbar links relative to page
        :param list[str] context: the context the links will be used in.
            The sequence of directories to be prepended to the HTML file in
            the resulting navbar
        :return str: navbar links as HTML-formatted string
        """
        # determine paths
        if wd is None and context is None:
            raise ValueError(
                "Either 'wd' (path the links should be relative to) or "
                "'context' (the context for the links) has to be provided."
            )
        status_relpath = _make_relpath(
            file_name=os.path.join(self.pipeline_reports, "status.html"),
            wd=wd,
            context=context,
        )
        objects_relpath = _make_relpath(
            file_name=os.path.join(self.pipeline_reports, "objects.html"),
            wd=wd,
            context=context,
        )
        samples_relpath = _make_relpath(
            file_name=os.path.join(self.pipeline_reports, "records.html"),
            wd=wd,
            context=context,
        )
        glossary_relpath = _make_relpath(
            file_name=os.path.join(self.pipeline_reports, "glossary.html"),
            wd=wd,
            context=context,
        )
        # determine the outputs IDs by type
        obj_result_ids = self.get_nonhighlighted_results(OBJECT_TYPES)

        # Remove project level objects because they are stored at the bottom of the index
        obj_to_remove = []
        for obj in obj_result_ids:
            if all(
                [
                    obj in list(self.prj.schema.project_level_data.keys()),
                    obj not in list(self.prj.schema.sample_level_data.keys()),
                ]
            ):
                obj_to_remove.append(obj)
        obj_result_ids = list(set(obj_to_remove) ^ set(obj_result_ids))

        dropdown_keys_objects = None
        dropdown_relpaths_objects = None
        sample_names = None
        if len(obj_result_ids) > 0:
            # If the number of objects is 20 or less, use a drop-down menu
            if len(obj_result_ids) <= 20:
                (
                    dropdown_relpaths_objects,
                    dropdown_keys_objects,
                ) = self._get_navbar_dropdown_data_objects(
                    objs=obj_result_ids, wd=wd, context=context
                )
        else:
            dropdown_relpaths_objects = objects_relpath
        if self.prj.record_count <= 20:
            (
                dropdown_relpaths_samples,
                sample_names,
            ) = self._get_navbar_dropdown_data_samples(wd=wd, context=context)
        else:
            # Create a menu link to the samples parent page
            dropdown_relpaths_samples = samples_relpath
        template_vars = dict(
            glossary_html_page=glossary_relpath,
            glossary_page_name="Glossary",
            status_html_page=status_relpath,
            status_page_name="Status",
            dropdown_keys_objects=dropdown_keys_objects,
            objects_page_name="Objects",
            samples_page_name="Records",
            objects_html_page=dropdown_relpaths_objects,
            samples_html_page=dropdown_relpaths_samples,
            menu_name_objects="Objects",
            menu_name_samples="Records",
            sample_names=sample_names,
            all_samples=samples_relpath,
            all_objects=objects_relpath,
            sample_reports_parent=None,
            project_report=project_index_html_relpath,
        )
        _LOGGER.debug(f"navbar_links.html | template_vars:\n{template_vars}")
        return render_jinja_template("navbar_links.html", self.jinja_env, template_vars)

    def create_object_htmls(self, navbar, footer):
        """
        Generates a page for an individual object type with all of its
        plots from each sample

        :param str navbar: HTML to be included as the navbar in the main summary page
        :param str footer: HTML to be included as the footer
        """
        file_results = self.get_nonhighlighted_results(["file"])
        image_results = self.get_nonhighlighted_results(["image"])

        if not os.path.exists(self.pipeline_reports):
            os.makedirs(self.pipeline_reports)
        for file_result in file_results:
            links = []
            html_page_path = os.path.join(
                self.pipeline_reports, f"{file_result}.html".replace(" ", "_").lower()
            )

            pipeline_types = ["sample", "project"]

            for pipeline_type in pipeline_types:
                self.prj.backend.pipeline_type = pipeline_type
                for sample in self.prj.backend.select_records()["records"]:
                    sample_name = sample["record_identifier"]
                    sample_result = fetch_pipeline_results(
                        project=self.prj,
                        sample_name=sample_name,
                    )
                    if file_result not in sample_result:
                        pass
                    else:
                        try:
                            if self.portable:
                                new_image_path = self._create_copy_for_porting(
                                    parent_path=sample_result[file_result]["path"],
                                    record_identifier=sample_name,
                                )
                                sample_result[file_result]["path"] = new_image_path

                            links.append(
                                [
                                    sample_name,
                                    os.path.relpath(
                                        sample_result[file_result]["path"],
                                        self.pipeline_reports,
                                    ),
                                ]
                            )
                        except Exception:
                            links.append(["LinkPathNotFound"])
                else:
                    link_desc = (
                        self.prj.result_schemas[file_result]["description"]
                        if "description" in self.prj.result_schemas[file_result]
                        else "No description in schema"
                    )
                    template_vars = dict(
                        navbar=navbar,
                        footer=footer,
                        name=file_result,
                        figures=[],
                        links=links,
                        desc=link_desc,
                        pipeline_name=self.pipeline_name,
                    )
                    save_html(
                        html_page_path,
                        render_jinja_template("object.html", self.jinja_env, args=template_vars),
                    )

        for image_result in image_results:
            html_page_path = os.path.join(self.pipeline_reports, f"{image_result}.html".lower())
            figures = []

            pipeline_types = ["sample", "project"]

            for pipeline_type in pipeline_types:
                self.prj.backend.pipeline_type = pipeline_type
                for sample in self.prj.backend.select_records()["records"]:
                    sample_name = sample["record_identifier"]
                    sample_result = fetch_pipeline_results(
                        project=self.prj,
                        sample_name=sample_name,
                    )
                    if image_result not in sample_result:
                        pass
                    else:
                        try:
                            if self.portable:
                                new_image_path = self._create_copy_for_porting(
                                    parent_path=sample_result[image_result]["path"],
                                    record_identifier=sample_name,
                                )

                                sample_result[image_result]["path"] = new_image_path

                                new_thumbnail_path = self._create_copy_for_porting(
                                    parent_path=sample_result[image_result]["thumbnail_path"],
                                    record_identifier=sample_name,
                                )

                                sample_result[image_result]["thumbnail_path"] = new_thumbnail_path

                            figures.append(
                                [
                                    os.path.relpath(
                                        sample_result[image_result]["path"],
                                        self.pipeline_reports,
                                    ),
                                    sample_name,
                                    os.path.relpath(
                                        sample_result[image_result]["thumbnail_path"],
                                        self.pipeline_reports,
                                    ),
                                ]
                            )
                        except:
                            figures.append(["FigurePathNotFound"])
                else:
                    img_desc = (
                        self.prj.result_schemas[image_result]["description"]
                        if "description" in self.prj.result_schemas[image_result]
                        else "No description in schema"
                    )
                    template_vars = dict(
                        navbar=navbar,
                        footer=footer,
                        name=image_result,
                        figures=figures,
                        links=[],
                        desc=img_desc,
                    )
                    _LOGGER.debug(f"object.html | template_vars:\n{template_vars}")
                    save_html(
                        html_page_path,
                        render_jinja_template("object.html", self.jinja_env, args=template_vars),
                    )
        self._reset_pipeline_type()

    def create_glossary_html(self, glossary_table, navbar, footer):
        template_vars = dict(
            glossary_table=glossary_table,
            navbar=navbar,
            footer=footer,
            pipeline_name=self.pipeline_name,
        )
        _LOGGER.debug(f"glossary.html | template_vars:\n{template_vars}")
        return render_jinja_template("glossary.html", self.jinja_env, template_vars)

    def create_sample_html(self, sample_stats, navbar, footer, sample_name):
        """
        Produce an HTML page containing all of a sample's objects
        and the sample summary statistics

        :param str sample_name: the name of the current sample
        :param dict sample_stats: pipeline run statistics for the current sample
        :param str navbar: HTML to be included as the navbar in the main summary page
        :param str footer: HTML to be included as the footer
        :param str pipeline_type: pipeline_type, 'project' or 'sample'
        :return str: path to the produced HTML page
        """
        if not os.path.exists(self.pipeline_reports):
            os.makedirs(self.pipeline_reports)
        html_page = os.path.join(
            self.pipeline_reports, f"{sample_name}.html".replace(" ", "_").lower()
        )
        if self.prj.cfg["multi_result_files"] is True:
            self.prj.cfg["record_identifier"] = sample_name
            temp_result_file_path = mkabs(
                self.prj.resolve_results_file_path(self.prj.cfg["unresolved_result_path"]),
                self.prj.cfg["config_path"],
            )
            make_subdirectories(temp_result_file_path)
            self.prj.backend.status_file_dir = os.path.dirname(
                mkabs(temp_result_file_path, self.prj.cfg["config_path"])
            )

        flag = self.prj.get_status(record_identifier=sample_name)

        if not flag:
            button_class = "btn btn-secondary"
            flag = "Missing"
        else:
            try:
                flag_dict = BUTTON_APPEARANCE_BY_FLAG[flag]
            except KeyError:
                button_class = "btn btn-secondary"
                flag = "Unknown"
            else:
                button_class = flag_dict["button_class"]
                flag = flag_dict["flag"]
        highlighted_results = fetch_pipeline_results(
            project=self.prj,
            sample_name=sample_name,
            inclusion_fun=lambda x: x == "file",
            highlighted=True,
        )

        for k in highlighted_results.keys():
            highlighted_results[k]["path"] = os.path.relpath(
                highlighted_results[k]["path"], self.pipeline_reports
            )

        links = []
        file_results = fetch_pipeline_results(
            project=self.prj,
            sample_name=sample_name,
            inclusion_fun=lambda x: x == "file",
        )
        for result_id, result in file_results.items():
            desc = (
                self.schema[result_id]["description"]
                if "description" in self.schema[result_id]
                else ""
            )
            if result:
                try:
                    if self.portable:
                        new_image_path = self._create_copy_for_porting(
                            parent_path=result["path"], record_identifier=sample_name
                        )
                        result["path"] = new_image_path

                    links.append(
                        [
                            f"<b>{result['title']}</b>: {desc}",
                            os.path.relpath(result["path"], self.pipeline_reports),
                        ]
                    )
                except FileNotFoundError:
                    _LOGGER.warning(f"File not found for {result_id} and {result}")
        image_results = fetch_pipeline_results(
            project=self.prj,
            sample_name=sample_name,
            inclusion_fun=lambda x: x == "image",
        )
        figures = []
        for result_id, result in image_results.items():
            if result:
                try:
                    if self.portable:
                        new_image_path = self._create_copy_for_porting(
                            parent_path=result["path"], record_identifier=sample_name
                        )

                        result["path"] = new_image_path

                        new_thumbnail_path = self._create_copy_for_porting(
                            parent_path=result["thumbnail_path"], record_identifier=sample_name
                        )

                        result["thumbnail_path"] = new_thumbnail_path

                    figures.append(
                        [
                            os.path.relpath(result["path"], self.pipeline_reports),
                            result["title"],
                            os.path.relpath(result["thumbnail_path"], self.pipeline_reports),
                        ]
                    )
                except FileNotFoundError:
                    _LOGGER.warning(f"File not found for {result_id} and {result}")

        template_vars = dict(
            report_class="Sample",
            navbar=navbar,
            footer=footer,
            sample_name=sample_name,
            links=links,
            figures=figures,
            button_class=button_class,
            sample_stats=sample_stats,
            flag=flag,
            highlighted_results=highlighted_results,
            pipeline_name=self.pipeline_name,
            amendments="",
        )
        _LOGGER.debug(f"sample.html | template_vars:\n{template_vars}")
        save_html(
            html_page,
            render_jinja_template("sample.html", self.jinja_env, template_vars),
        )
        return html_page

    def create_status_html(self, status_table, navbar, footer):
        """
        Generates a page listing all the samples, their run status, their
        log file, and the total runtime if completed.

        :param str navbar: HTML to be included as the navbar in the main summary page
        :param str footer: HTML to be included as the footer
        :return str: rendered status HTML file
        """
        _LOGGER.debug("Building status page...")
        template_vars = dict(
            status_table=status_table,
            navbar=navbar,
            footer=footer,
            pipeline_name=self.pipeline_name,
        )
        _LOGGER.debug(f"status.html | template_vars:\n{template_vars}")
        return render_jinja_template("status.html", self.jinja_env, template_vars)

    def create_index_html(self, navbar, footer):
        """
        Generate an index.html style project home page w/ sample summary
        statistics

        :param str navbar: HTML to be included as the navbar in the main
            summary page
        :param str footer: HTML to be included as the footer
        """
        # set default encoding when running in python2
        if sys.version[0] == "2":
            from importlib import reload

            reload(sys)
            sys.setdefaultencoding("utf-8")
        _LOGGER.info(f"Building index page for pipeline: {self.pipeline_name}")

        # Create stats and object summaries
        table_path_list = _create_stats_objs_summaries(self.prj, self.pipeline_name)

        # Add stats_summary.tsv button link
        stats_file_path = table_path_list[0]

        stats_file_path = stats_file_path if os.path.exists(stats_file_path) else None

        # Add objects_summary.yaml button link

        objs_file_path = table_path_list[1]

        objs_file_path = objs_file_path if os.path.exists(objs_file_path) else None

        if self.portable:
            stats_file_path = self._create_copy_for_porting(stats_file_path, "")
            stats_file_path = os.path.relpath(
                stats_file_path,
                self.pipeline_reports,
            )
            objs_file_path = self._create_copy_for_porting(objs_file_path, "")
            objs_file_path = os.path.relpath(
                objs_file_path,
                self.pipeline_reports,
            )

        # Add stats summary table to index page and produce individual
        # sample pages
        # Produce table rows
        table_row_data = []
        _LOGGER.info(" * Creating sample pages")

        # get all potentially reportable keys from the output schema
        all_result_identifiers = list(self.schema.keys())

        if self.looper_samples is not None:
            # If looper passes the samples from the PEP, we should add the attributes to the cell tables:
            input_sample_attributes = self.looper_samples[0]._mapped_attr["_attributes"]
            # Place at front of columns
            all_result_identifiers = input_sample_attributes + all_result_identifiers

        if self.prj.cfg["multi_result_files"] is True:
            pipeline_types = ["sample"]
        else:
            pipeline_types = [self.prj.backend.pipeline_type]

        sorted_sample_stat_results = {}
        for pipeline_type in pipeline_types:
            self.prj.backend.pipeline_type = pipeline_type
            for sample in self.prj.backend.select_records()["records"]:
                sample_name = sample["record_identifier"]
                sample_stat_results = fetch_pipeline_results(
                    project=self.prj,
                    sample_name=sample_name,
                    inclusion_fun=None,
                    casting_fun=str,
                )

                if self.looper_samples is not None:
                    for s in range(len(self.looper_samples)):
                        if sample_name == self.looper_samples[s]._mapped_attr["sample_name"]:
                            for attribute in input_sample_attributes:
                                value = self.looper_samples[s]._mapped_attr[attribute]
                                sample_stat_results.update({attribute: value})

                for key in all_result_identifiers:
                    if key not in sample_stat_results.keys():
                        sample_stat_results[key] = ""

                for key in self.schema.keys():
                    if "type" in self.schema[key]:
                        if (
                            self.schema[key]["type"] == "file"
                            or self.schema[key]["type"] == "image"
                            or self.schema[key]["type"] == "object"
                        ):
                            del sample_stat_results[key]

                # Sort to ensure alignment in the table
                sorted_sample_stat_results = dict(sorted(sample_stat_results.items()))

                sample_html = self.create_sample_html(
                    sorted_sample_stat_results,
                    navbar,
                    footer,
                    sample_name,
                )
                rel_sample_html = os.path.relpath(sample_html, self.pipeline_reports)
                # treat sample_name column differently - will need to provide
                # a link to the sample page
                table_cell_data = [[rel_sample_html, sample_name]]
                table_cell_data += list(sorted_sample_stat_results.values())
                table_row_data.append(table_cell_data)
        self._reset_pipeline_type()
        # Create parent samples page with links to each sample
        save_html(
            path=os.path.join(self.pipeline_reports, "records.html"),
            template=self.create_sample_parent_html(navbar, footer),
        )
        _LOGGER.info(" * Creating object pages")
        # Create objects pages
        self.create_object_htmls(navbar, footer)

        # Create parent objects page with links to each object type
        save_html(
            path=os.path.join(self.pipeline_reports, "objects.html"),
            template=self.create_object_parent_html(navbar, footer),
        )
        # Create status page with each sample's status listed
        status_tab = create_status_table(
            report_obj=self,
            project=self.prj,
            pipeline_reports_dir=self.pipeline_reports,
            portable=self.portable,
        )
        save_html(
            path=os.path.join(self.pipeline_reports, "status.html"),
            template=self.create_status_html(status_tab, navbar, footer),
        )
        # glossary_html = self.create_glossary_html(navbar,footer)
        glossary_table = create_glossary_table(project=self.prj)
        save_html(
            path=os.path.join(self.pipeline_reports, "glossary.html"),
            template=self.create_glossary_html(glossary_table, navbar, footer),
        )

        project_objects = self.create_project_objects()
        columns_table = ["Record Identifiers"] + list(sorted_sample_stat_results.keys())
        columns_stats = list(sorted_sample_stat_results.keys())
        template_vars = dict(
            navbar=navbar,
            stats_file_path=stats_file_path,
            objs_file_path=objs_file_path,
            columns=columns_table,
            columns_json=dumps(columns_stats),
            table_row_data=table_row_data,
            project_name=self.prj.cfg[PROJECT_NAME],
            pipeline_name=self.pipeline_name,
            stats_json=self._stats_to_json_str(),
            project_objects=project_objects,
            footer=footer,
            amendments="",
        )
        _LOGGER.debug(f"index.html | template_vars:\n{template_vars}")
        save_html(
            self.index_html_path,
            render_jinja_template("index.html", self.jinja_env, template_vars),
        )

    def create_project_objects(self):
        """
        Render available project level outputs defined in the
        pipeline output schemas
        """
        _LOGGER.debug("Building project objects section...")
        figures = []
        links = []
        warnings = []

        file_results = self.get_nonhighlighted_results(["file"])
        image_results = self.get_nonhighlighted_results(["image"])

        if not os.path.exists(self.pipeline_reports):
            os.makedirs(self.pipeline_reports)
        for file_result in file_results:
            html_page_path = os.path.join(self.pipeline_reports, f"{file_result}.html".lower())

            pipeline_types = ["project"]

            for pipeline_type in pipeline_types:
                self.prj.backend.pipeline_type = pipeline_type
                for sample in self.prj.backend.select_records()["records"]:
                    sample_name = sample["record_identifier"]
                    sample_result = fetch_pipeline_results(
                        project=self.prj,
                        sample_name=sample_name,
                    )
                    if file_result not in sample_result:
                        pass
                    else:
                        try:
                            if self.portable:
                                new_image_path = self._create_copy_for_porting(
                                    parent_path=sample_result[file_result]["path"],
                                    record_identifier=sample_name,
                                )
                                sample_result[file_result]["path"] = new_image_path

                            links.append(
                                [
                                    f"<b>{sample_result[file_result]['title']}</b>",
                                    os.path.relpath(
                                        sample_result[file_result]["path"],
                                        self.pipeline_reports,
                                    ),
                                ]
                            )
                        except Exception:
                            links.append(["LinkPathNotFound"])
                else:
                    link_desc = (
                        self.prj.result_schemas[file_result]["description"]
                        if "description" in self.prj.result_schemas[file_result]
                        else "No description in schema"
                    )

        for image_result in image_results:
            html_page_path = os.path.join(self.pipeline_reports, f"{image_result}.html".lower())

            pipeline_types = ["project"]

            for pipeline_type in pipeline_types:
                self.prj.backend.pipeline_type = pipeline_type
                for sample in self.prj.backend.select_records()["records"]:
                    sample_name = sample["record_identifier"]
                    sample_result = fetch_pipeline_results(
                        project=self.prj,
                        sample_name=sample_name,
                    )
                    if image_result not in sample_result:
                        pass
                    else:
                        try:
                            if self.portable:
                                new_image_path = self._create_copy_for_porting(
                                    parent_path=sample_result[image_result]["path"],
                                    record_identifier=sample_name,
                                )

                                sample_result[image_result]["path"] = new_image_path

                                new_thumbnail_path = self._create_copy_for_porting(
                                    parent_path=sample_result[image_result]["thumbnail_path"],
                                    record_identifier=sample_name,
                                )

                                sample_result[image_result]["thumbnail_path"] = new_thumbnail_path

                            figures.append(
                                [
                                    os.path.relpath(
                                        sample_result[image_result]["path"],
                                        self.pipeline_reports,
                                    ),
                                    sample_result[image_result]["title"],
                                    os.path.relpath(
                                        sample_result[image_result]["thumbnail_path"],
                                        self.pipeline_reports,
                                    ),
                                ]
                            )
                        except:
                            figures.append(["FigurePathNotFound"])
                else:
                    img_desc = (
                        self.prj.result_schemas[image_result]["description"]
                        if "description" in self.prj.result_schemas[image_result]
                        else "No description in schema"
                    )
        self._reset_pipeline_type()
        template_vars = dict(figures=figures, links=links)
        return render_jinja_template("project_object.html", self.jinja_env, template_vars)

    def get_nonhighlighted_results(self, types):
        """
        Get a list of non-highlighted results in the schema

        :param list[str] types: types to narrow down the results
        :return list[str]: result ID that are of the requested type and
            are not highlighted
        """
        results = []

        for key, value in self.schema.items():
            if self.schema[key]["type"] in types:
                if "highlight" not in self.schema[key].keys():
                    results.append(key)
                # intentionally "== False" to exclude "falsy" values
                elif self.schema[key]["highlight"] is False:
                    results.append(key)

        return results

    def _stats_to_json_str(self):
        results = {}
        if self.prj.cfg["multi_result_files"] is True:
            pipeline_types = ["sample"]
        else:
            pipeline_types = [self.prj.backend.pipeline_type]

        for pipeline_type in pipeline_types:
            self.prj.backend.pipeline_type = pipeline_type
            for sample in self.prj.backend.select_records()["records"]:
                sample_name = sample["record_identifier"]
                results[sample_name] = fetch_pipeline_results(
                    project=self.prj,
                    sample_name=sample_name,
                    inclusion_fun=lambda x: x not in OBJECT_TYPES,
                    casting_fun=str,
                )
        self._reset_pipeline_type()
        return dumps(results)

    def _get_navbar_dropdown_data_objects(self, objs, wd, context):
        if objs is None or len(objs) == 0:
            return None, None
        relpaths = []
        displayable_ids = []
        for obj_id in objs:
            displayable_ids.append(obj_id.replace("_", " "))
            page_name = os.path.join(
                self.pipeline_reports, (obj_id + ".html").replace(" ", "_").lower()
            )
            relpaths.append(_make_relpath(page_name, wd, context))
        return relpaths, displayable_ids

    def _get_navbar_dropdown_data_samples(self, wd, context):
        relpaths = []
        sample_names = []
        if self.prj.cfg["multi_result_files"] is True:
            pipeline_types = ["sample", "project"]
        else:
            pipeline_types = [self.prj.backend.pipeline_type]

        for pipeline_type in pipeline_types:
            self.prj.backend.pipeline_type = pipeline_type
            for sample in self.prj.backend.select_records()["records"]:
                sample_name = sample["record_identifier"]
                page_name = os.path.join(
                    self.pipeline_reports,
                    f"{sample_name}.html".replace(" ", "_").lower(),
                )
                relpaths.append(_make_relpath(page_name, wd, context))
                sample_names.append(sample_name)
        self._reset_pipeline_type()
        return relpaths, sample_names


def render_jinja_template(name, jinja_env, args=dict()):
    """
    Render template in the specified jinja environment using the provided args

    :param str name: name of the template
    :param dict args: arguments to pass to the template
    :param jinja2.Environment jinja_env: the initialized environment to use in
        this the looper HTML reports context
    :return str: rendered template
    """
    assert isinstance(args, dict), "args has to be a dict"
    template = jinja_env.get_template(name)
    return template.render(**args)


def save_html(path, template):
    """
    Save rendered template as an HTML file

    :param str path: the desired location for the file to be produced
    :param str template: the template or just string
    """
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    try:
        with open(path, "w") as f:
            f.write(template)
    except IOError:
        _LOGGER.error("Could not write the HTML file: {}".format(path))


def get_jinja_env(templates_dirname=None):
    """
    Create jinja environment with the provided path to the templates directory

    :param str templates_dirname: path to the templates directory
    :return jinja2.Environment: jinja environment
    """
    if templates_dirname is None:
        file_dir = os.path.dirname(os.path.realpath(__file__))
        templates_dirname = os.path.join(file_dir, TEMPLATES_DIRNAME)
    _LOGGER.debug("Using templates dir: " + templates_dirname)
    return jinja2.Environment(loader=jinja2.FileSystemLoader(templates_dirname))


def _get_file_for_sample(prj, sample_name, appendix, pipeline_name=None, basename=False):
    """
    Safely looks for files matching the appendix in the specified
    location for the sample

    :param str sample_name: name of the sample that the file name
        should be found for
    :param str appendix: the ending pecific for the file
    :param bool basename: whether to return basename only
    :return str: the name of the matched file
    """
    fp = os.path.join(prj.results_folder, sample_name)
    prepend_name = ""
    if pipeline_name:
        prepend_name += pipeline_name
    if hasattr(prj, AMENDMENTS_KEY) and getattr(prj, AMENDMENTS_KEY):
        prepend_name += f"_{'_'.join(getattr(prj, AMENDMENTS_KEY))}"
    prepend_name = prepend_name + "_" if prepend_name else ""
    fp = os.path.join(fp, f"{prepend_name}{appendix}")
    if os.path.exists(fp):
        return os.path.basename(fp) if basename else fp
    raise FileNotFoundError(fp)


def _get_relpath_to_file(file_name, sample_name, location, relative_to):
    """
    Safely gets the relative path for the file for the specified sample

    :param str file_name: name of the file
    :param str sample_name: name of the sample that the file path
        should be found for
    :param str location: where to look for the file
    :param str relative_to: path the result path should be relative to
    :return str: a path to the file
    """
    abs_file_path = os.path.join(location, sample_name, file_name)
    rel_file_path = os.path.relpath(abs_file_path, relative_to)
    if file_name is None or not os.path.exists(abs_file_path):
        return None
    return rel_file_path


def _make_relpath(file_name, wd, context=None):
    """
    Create a path relative to the context. This function introduces the
    flexibility to the navbar links creation, which the can be used outside
    of the native looper summary pages.

    :param str file_name: the path to make relative
    :param str wd: the dir the path should be relative to
    :param list[str] context: the context the links will be used in. The
        sequence of directories to be prepended to the HTML
        file in the resulting navbar
    :return str: relative path
    """
    relpath = os.path.relpath(file_name, wd)
    return relpath if not context else os.path.join(os.path.join(*context), relpath)


def _read_csv_encodings(path, encodings=["utf-8", "ascii"], **kwargs):
    """
    Try to read file with the provided encodings

    :param str path: path to file
    :param list encodings: list of encodings to try
    """
    idx = 0
    while idx < len(encodings):
        e = encodings[idx]
        try:
            t = _pd.read_csv(path, encoding=e, **kwargs)
            return t
        except UnicodeDecodeError:
            pass
        idx = idx + 1
    _LOGGER.warning(f"Could not read the log file '{path}' with encodings '{encodings}'")


def _read_tsv_to_json(path):
    """
    Read a tsv file to a JSON formatted string

    :param path: to file path
    :return str: JSON formatted string
    """
    assert os.path.exists(path), "The file '{}' does not exist".format(path)
    _LOGGER.debug("Reading TSV from '{}'".format(path))
    df = _pd.read_csv(path, sep="\t", index_col=False, header=None)
    return df.to_json()


def fetch_pipeline_results(
    project,
    sample_name=None,
    inclusion_fun=None,
    casting_fun=None,
    highlighted=False,
):
    """
    Get the specific pipeline results for sample based on inclusion function

    :param looper.Project project: project to get the results for
    :param str sample_name: sample ID
    :param callable(str) inclusion_fun: a function that determines whether the
        result should be returned based on it's type. Example input that the
        function will be fed with is: 'image' or 'integer'
    :param callable(str) casting_fun: a function that will be used to cast the
        each of the results to a proper type before returning, e.g int, str
    :param bool highlighted: return the highlighted or regular results
    :param str pipeline_type: pipeline_type, 'project' or 'sample'
    :return dict: selected pipeline results
    """

    def pass_all_fun(x):
        return x

    inclusion_fun = inclusion_fun or pass_all_fun
    casting_fun = casting_fun or pass_all_fun
    psm = project
    # exclude object-like results from the stats results mapping
    rep_data = psm.retrieve_one(record_identifier=sample_name)
    results = {
        k: casting_fun(v)
        for k, v in rep_data.items()
        if k in psm.result_schemas
        and inclusion_fun(psm.result_schemas[k].get("object_type", psm.result_schemas[k]["type"]))
    }
    if highlighted:
        return {k: v for k, v in results.items() if k in psm.highlighted_results}
    return {k: v for k, v in results.items() if k not in psm.highlighted_results}


def uniqify(seq):
    """Fast way to uniqify while preserving input order."""
    # http://stackoverflow.com/questions/480214/
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def create_status_table(report_obj, project, pipeline_reports_dir: str, portable: bool) -> str:
    """
    Creates status table, the core of the status page.

    :param report_obj: The HTML builder object
    :param PipestatManager project: project to get the results for
    :param str pipeline_reports_dir: path to the pipeline reports directory
    :param bool portable: is the report to be portable?
    :return str: rendered status HTML file
    """

    def _rgb2hex(r, g, b):
        return "#{:02x}{:02x}{:02x}".format(r, g, b)

    def _warn(what, e, sn):
        _LOGGER.debug(
            f"Caught exception: {e}\n"
            f"Could not determine {what} for sample: {sn}. "
            f"Not reported or pipestat status schema is faulty."
        )

    log_paths = []
    log_link_names = []
    sample_paths = []
    sample_names = []
    statuses = []
    status_styles = []
    times = []
    mems = []
    status_descs = []

    if project.cfg["multi_result_files"] is True:
        pipeline_types = ["sample", "project"]
    else:
        pipeline_types = [project.backend.pipeline_type]

    for pipeline_type in pipeline_types:
        project.backend.pipeline_type = pipeline_type
        for sample in project.backend.select_records()["records"]:
            sample_name = sample["record_identifier"]
            psm = project
            sample_names.append(sample_name)
            # status and status style
            try:
                if psm.cfg["multi_result_files"] is True:
                    psm.cfg["record_identifier"] = sample_name
                    temp_result_file_path = mkabs(
                        psm.resolve_results_file_path(psm.cfg["unresolved_result_path"]),
                        psm.cfg["config_path"],
                    )
                    make_subdirectories(temp_result_file_path)
                    psm.backend.status_file_dir = os.path.dirname(
                        mkabs(temp_result_file_path, psm.cfg["config_path"])
                    )
                status = psm.get_status(record_identifier=sample_name)
                statuses.append(status)
                status_metadata = psm.status_schema[status]
                status_styles.append(_rgb2hex(*status_metadata["color"]))
                status_descs.append(status_metadata["description"])
            except Exception as e:
                _warn("status", e, sample_name)
                statuses.append(NO_DATA_PLACEHOLDER)
                status_styles.append(NO_DATA_PLACEHOLDER)
                status_descs.append(NO_DATA_PLACEHOLDER)
            sample_paths.append(f"{sample_name}.html".replace(" ", "_").lower())
            # log file path
            try:
                log = glob.glob(psm.backend.status_file_dir + "**/*log.md")[
                    0
                ]  # Assumes the log file will be in status dir
                assert os.path.exists(log), FileNotFoundError(f"Not found: {log}")
                if portable:
                    new_log_path = report_obj._create_copy_for_porting(
                        parent_path=log,
                        record_identifier=sample_name,
                    )
                    log = new_log_path
                log_link_names.append(os.path.basename(log))
                log_paths.append(os.path.relpath(log, pipeline_reports_dir))
            except Exception as e:
                _warn("log", e, sample)
                log_link_names.append(NO_DATA_PLACEHOLDER)
                log_paths.append("")
            # runtime and peak mem
            try:
                profile = glob.glob(psm.backend.status_file_dir + "**/*profile.tsv")[
                    0
                ]  # Assumes the profile file will be in status dir
                assert os.path.exists(profile), FileNotFoundError(f"Not found: {profile}")
                df = _pd.read_csv(profile, sep="\t", comment="#", names=PROFILE_COLNAMES)
                df["runtime"] = _pd.to_timedelta(df["runtime"])
                times.append(_get_runtime(df))
                mems.append(_get_maxmem(df))
            except Exception as e:
                _warn("profile", e, sample)
                times.append(NO_DATA_PLACEHOLDER)
                mems.append(NO_DATA_PLACEHOLDER)

    project.backend.pipeline_type = project.pipeline_type

    template_vars = dict(
        sample_names=sample_names,
        log_paths=log_paths,
        status_styles=status_styles,
        statuses=statuses,
        times=times,
        mems=mems,
        sample_paths=sample_paths,
        log_link_names=log_link_names,
        status_descs=status_descs,
    )
    _LOGGER.debug(f"status_table.html | template_vars:\n{template_vars}")
    return render_jinja_template("status_table.html", get_jinja_env(), template_vars)


def create_glossary_table(project):
    items = []
    descs = []

    for k, v in project.result_schemas.items():
        items.append(k)
        descs.append(v["description"])

    template_vars = dict(
        items=items,
        descriptions=descs,
    )

    # return (items,descs)
    return render_jinja_template("glossary_table.html", get_jinja_env(), template_vars)


def _get_maxmem(profile: _pd.DataFrame) -> str:
    """
    Get current peak memory

    :param pandas.DataFrame profile: a data frame representing
        the current profile.tsv for a sample
    :return str: max memory
    """
    return f"{str(max(profile['mem']) if not profile['mem'].empty else 0)} GB"


def _get_runtime(profile_df: _pd.DataFrame) -> str:
    """
    Collect the unique and last duplicated runtimes, sum them and then
    return in str format

    :param pandas.DataFrame profile_df: a data frame representing
        the current profile.tsv for a sample
    :return str: sum of runtimes
    """
    unique_df = profile_df[~profile_df.duplicated("cid", keep="last").values]
    return str(
        timedelta(seconds=sum(unique_df["runtime"].apply(lambda x: x.total_seconds())))
    ).split(".")[0]


def get_file_for_table(prj, pipeline_name: str, appendix=None, directory=None) -> str:
    """
    Create a path to the file for the current project.
    Takes the possibility of amendment being activated at the time

    Format of the output path:
    {output_dir}/{directory}/{p.name}_{pipeline_name}_{active_amendments}_{appendix}

    :param pipestat manager object prj: project object
    :param str pipeline_name: name of the pipeline to get the file for
    :param str appendix: the appendix of the file to create the path for,
        like 'objs_summary.tsv' for objects summary file
    :param directory: subdirectory (if desired)
    :return str fp: path to the file
    """
    # TODO make determining the output_dir its own small function since we use the same code in HTML report building.
    results_file_path = getattr(prj.backend, "results_file_path", None)
    config_path = prj.cfg.get("config_path", None)
    output_dir = prj.cfg.get(OUTPUT_DIR, None)
    table_dir = output_dir or results_file_path or config_path
    if not os.path.isdir(table_dir):
        table_dir = os.path.dirname(table_dir)
    fp = os.path.join(table_dir, directory or "", f"{pipeline_name}")
    if hasattr(prj, "amendments") and getattr(prj, "amendments"):
        fp += f"_{'_'.join(prj.amendments)}"
    if appendix:
        fp += f"_{appendix}"
    return fp


def _create_stats_objs_summaries(prj, pipeline_name: str) -> List[str]:
    """
    Create stats spreadsheet and objects summary.

    :param pipestat.PipestatManager prj: pipestat object used to create table
    :param str pipeline_name: name of the pipeline to tabulate results for
    :return List[str] [tsv_outfile_path, objs_yaml_path]: list of paths to tsv_outfile_path, objs_yaml_path
    """

    _LOGGER.info("Creating objects summary")
    reported_objects = {}
    stats = []

    columns = ["Sample Index", "Record Identifier", "Pipeline Type"]

    record_index = 0

    all_result_identifiers = list(prj.result_schemas.keys())

    if prj.cfg["multi_result_files"] is True:
        pipeline_types = ["sample", "project"]
    else:
        pipeline_types = [prj.backend.pipeline_type]

    for pipeline_type in pipeline_types:
        prj.backend.pipeline_type = pipeline_type
        records = prj.backend.select_records()["records"]
        for record in records:
            record_index += 1
            record_name = record["record_identifier"]

            if prj.cfg[PIPELINE_TYPE] == "sample":
                reported_stats = [record_index, record_name, "sample"]
                rep_data = prj.retrieve_one(record_identifier=record_name)
            else:
                rep_data = prj.retrieve_one(record_identifier=record_name)
                reported_stats = [record_index, record_name, "project"]
            for key in all_result_identifiers:
                if key not in rep_data.keys():
                    rep_data[key] = "None reported"

            # Sort to ensure alignment in the table
            rep_data = dict(sorted(rep_data.items()))

            for k, v in rep_data.items():
                if v:
                    if k in prj.result_schemas and prj.result_schemas[k]["type"] in OBJECT_TYPES:
                        if k in all_result_identifiers:
                            all_result_identifiers.remove(k)
                        if v != "None reported":
                            if isinstance(v, list):
                                sample_reported_objects = {k: v}
                            else:
                                sample_reported_objects = {k: dict(v)}
                        else:
                            sample_reported_objects = {k: "None reported"}
                        if record_name in reported_objects:
                            reported_objects[record_name].update(sample_reported_objects)
                        else:
                            reported_objects.update({record_name: sample_reported_objects})
                    if (
                        k in prj.result_schemas
                        and prj.result_schemas[k]["type"] not in OBJECT_TYPES
                    ):
                        if k not in columns:
                            columns.append(k)
                        reported_stats.append(v)

            stats.append(reported_stats)
    prj.backend.pipeline_type = prj.pipeline_type
    # Stats File
    tsv_outfile_path = get_file_for_table(prj, pipeline_name, "stats_summary.tsv")
    stats.insert(0, columns)
    with open(tsv_outfile_path, "w") as tsv_outfile:
        writer = csv.writer(tsv_outfile, delimiter="\t")
        for row in stats:
            writer.writerow(row)
    _LOGGER.info(f"'{pipeline_name}' pipeline stats summary (n={len(stats)}): {tsv_outfile_path}")

    # Create Object yaml
    objs_yaml_path = get_file_for_table(prj, pipeline_name, "objs_summary.yaml")
    with open(objs_yaml_path, "w") as outfile:
        yaml.dump(reported_objects, outfile)
    _LOGGER.info(
        f"'{pipeline_name}' pipeline objects summary (n={len(reported_objects.keys())}): {objs_yaml_path}"
    )

    return [tsv_outfile_path, objs_yaml_path]
