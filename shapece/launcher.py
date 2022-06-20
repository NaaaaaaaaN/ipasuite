#! /usr/bin/env python3

import snakemake as sm
import os
import shutil
import fire
import subprocess

base_path = os.path.dirname(__file__)


class Launcher(object):
    def __init__(self, config: str = None, cores: int = 8, stoponerror: bool = False):
        if config is not None:
            config = [config]
        self._config = config
        self._cores = cores
        self._keepgoing = not stoponerror

    def config(self):
        path = os.path.join(base_path, "configurator.ipynb")
        subprocess.run(["voila", path])

    def init(self, project: str):
        if os.path.exists(project):
            fire.core.FireError(f"{project} folder already exists")

        os.makedirs(os.path.join(project, "resources"))
        os.mkdir(os.path.join(project, "results"))
        os.mkdir(os.path.join(project, "config"))
        shutil.copy(
            os.path.join(base_path, "config", "config.tpl.yaml"),
            os.path.join(project, "config", "config.yaml"),
        )
        shutil.copy(
            os.path.join(base_path, "config", "samples.tpl.tsv"),
            os.path.join(project, "config", "samples.tsv"),
        )

    def refactor(self, action: str = "addpositions"):
        extra_config = dict()
        if action == "addpositions":
            extra_config["refactor_enabled"] = True
            targets = ["all_add_positions"]
        else:
            fire.core.FireError(f"invalid refactor option {action}")

        try:
            sm.snakemake(
                os.path.join(base_path, "workflow", "Snakefile"),
                configfiles=self._config,
                config=extra_config,
                targets=targets,
                cores=self._cores,
                keepgoing=self._keepgoing,
                use_conda=True,
                conda_prefix="~/.shapece/conda",
            )
        except Exception as e:
            print(e)

    def run(
        self,
        action="all",
    ):

        targets = ["all"]
        extra_config = dict()

        try:
            sm.snakemake(
                os.path.join(base_path, "workflow", "Snakefile"),
                configfiles=self._config,
                config=extra_config,
                targets=targets,
                cores=self._cores,
                keepgoing=self._keepgoing,
                use_conda=True,
                conda_prefix="~/.shapece/conda",
            )
        except Exception as e:
            print(e)


def main_wrapper():
    fire.Fire(Launcher)
