#! /usr/bin/env python3

import snakemake as sm
from snakemake.utils import validate
import os
import shutil
import fire
import subprocess
from glob import glob
from ruamel.yaml import YAML
import logging
from .workflow.rules import load_samples
from .workflow.scripts.tools.utils import fasta
import filecmp
import pathlib
from colorlog import ColoredFormatter
import pandas as pd
import traceback

import csv
import re
import RNA
import varnaapi

LOGFORMAT = "  %(log_color)s%(levelname)-8s%(reset)s | %(message)s"
DATE_FORMAT = "%b %d %H:%M:%S"
formatter = ColoredFormatter(LOGFORMAT)
handler = logging.StreamHandler()
handler.setFormatter(formatter)
# logging.basicConfig(format="%(levelname)s:%(message)s")
logger = logging.getLogger(__name__)
logger.addHandler(handler)
yaml = YAML()
base_path = os.path.dirname(__file__)

step_order = [
    "fluo-ceq8000",
    "fluo-ce",
    "subseq",
    "qushape",
    "reactivity",
    "normreact",
    "alignnormreact",
    "aggreact",
    "aggreact-ipanemap",
    "ipanemap-config",
    "ipanemap-out",
    "structure",
    "varna",
    "footprint",
]


class Launcher(object):
    def __init__(
        self,
        config: str = None,
        cores: int = 8,
        stoponerror: bool = False,
        verbose: bool = False,
    ):
        self._config = config
        self._cores = cores
        self._keepgoing = not stoponerror
        self._verbose = verbose

    def _choose_config(self, config):
        config = (
            config
            if config is not None
            else (
                "config.yaml" if os.path.exists("config.yaml") else "config/config.yaml"
            )
        )
        if not os.path.exists(config):
            raise fire.core.FireError(
                f"{config} file does not exist, please init your "
                "project using `ipasuite init` command or specify another path"
            )
        return config

    def config(self, dev=False):
        """
        Configure the pipeline using the configurator.
        """
        self._config = self._choose_config(self._config)
        path = os.path.join(base_path, "configurator.ipynb")
        env = os.environ.copy()
        env["CONFIG_FILE_PATH"] = os.path.join(self._config)
        env["PROJECT_PATH"] = os.path.join(os.getcwd())
        if dev:
            subprocess.Popen(["jupyter-notebook", path], env=env)
        else:
            subprocess.Popen(["voila", path], env=env)

    #    def report(self):
    #        targets = ["all"]
    #        extra_config = dict()
    #        report_path =os.path.join(os.getcwd(),"report.html")
    #        try:
    #            sm.snakemake(
    #                os.path.join(base_path, "workflow", "Snakefile"),
    #                configfiles=[self._config],
    #                config=extra_config,
    #                targets=targets,
    #                cores=self._cores,
    #                report=report_path,
    #                keepgoing=self._keepgoing,
    #                use_conda=True,
    #                verbose=self._verbose,
    #                conda_prefix="~/.ipasuite/conda",
    #            )
    #        except Exception as e:
    #            print(e)
    #        webbrowser.open_new_tab(report_path)
    def report(self, dev=False):
        """
        Generate a summary report of results.
        """
        self._config = self._choose_config(self._config)
        path = os.path.join(base_path, "report.ipynb")
        env = os.environ.copy()
        env["CONFIG_FILE_PATH"] = os.path.join(self._config)
        env["PROJECT_PATH"] = os.path.join(os.getcwd())
        if dev:
            subprocess.Popen(["jupyter-notebook", path], env=env)
        else:
            subprocess.Popen(["voila", path], env=env)

    def init(self, project: str):
        """
        Generate the project folder.
        """
        if os.path.exists(project):
            fire.core.FireError(f"{project} folder already exists")
        os.mkdir(project)
        shutil.copy(
            os.path.join(base_path, "config", "config.tpl.yaml"),
            os.path.join(project, "config.yaml"),
        )
        shutil.copy(
            os.path.join(base_path, "config", "samples.tpl.tsv"),
            os.path.join(project, "samples.tsv"),
        )

        os.makedirs(os.path.join(project, "resources/raw_data"))
        os.mkdir(os.path.join(project, "results"))


    def prep(self):
        """
        Autofill samples.tsv file.
        """
        response = input("Do you want to overwrite sample.tsv? Please enter yes or no:")
        while response.lower() not in {'yes', 'no'}:
            print("Please enter yes or no")
            response = input("Do you want to overwrite sample.tsv? Please enter yes or no:")
            
        if response.lower() == 'yes':
            sample_file = os.path.join('samples.tsv')

            # copy existing configuration to backup file
            mydir = os.listdir('./')
            indices = []
            for file in mydir:
                m = re.match(rf'{sample_file}.(\d+)',file)
                if m:
                    indices.append(int(m[1]))
            bkp_file_index = max([0]+indices)+1
            bkp_file = f"{sample_file}.{bkp_file_index}"
            shutil.copy(sample_file, bkp_file)

            data_folder = os.path.join("resources/raw_data")
            files = os.listdir(data_folder)

            with open(sample_file, 'r', newline='') as sample:
                reader = csv.reader(sample, delimiter='\t')
                rows = next(reader)

            data = pd.DataFrame(columns=['rna_id', 'date', 'experimenter','probe', 'temperature', 'magnesium', 'ddNTP', 'rna_begin', 'rna_end',\
                'rt_begin_pos', 'rt_end_pos', 'replicate', 'probe_file', 'control_file', 'qushape_file','reference_qushape_file','map_file','discard'])

            template = ['{rna_id}_{probe_control_flag_unique_experiment_id}_{probe}_{temperature}_{magnesium}_{ddNTP}_{date}_{experimenter}.*.txt',\
                '{rna_id}_{probe}_{temperature}_{magnesium}_{replicate}.*.qushapey','{rna_id}_{probe}_{temperature}_{magnesium}_{replicate}.*.map']

            def translate_template_to_re(template):
                template = re.sub(r'\{([^}]+)\}', r'(?P<\1>[^_\.]+)',template)
                return template

            patterns = [translate_template_to_re(t) for t in template]

            file_info = []
            for file in files:
                row = {}
                m=None
                for p in patterns:
                    m = re.match(p,file)
                    if m is not None: break
                if m is not None:
                    list_info = m.groupdict()
                    for info in list(data):
                        if info in list_info:
                            row[info] = list_info[info]
                        else:
                            row[info] = ''
                    row['temperature'] = list_info['temperature'][1:]
                    if file.split('.')[-1] == 'txt' and list_info['probe_control_flag_unique_experiment_id'][0] == 'B':
                        replicate = file.split('_')[1][1:]
                        control_file = [file_name for file_name in files if len(file_name.split('_')) >1 and replicate == file_name.split('_')[1][1:] and \
                                        list_info['magnesium'] == file_name.split('_')[4] and list_info['temperature'] == file_name.split('_')[3] and \
                                        len(file_name.split('_')[2].split('-')) >1 and list_info['probe'] == file_name.split('_')[2].split('-')[1] and file_name != file][0]
                        row['probe_file'] = file
                        row['control_file'] = control_file
                        row['replicate'] = replicate
                        file_info.append(row)
                    elif file.split('.')[-1] == 'qushapey':
                        row['qushape_file'] = file
                        file_info.append(row)
                    elif file.split('.')[-1] == 'map':
                        row['map_file'] = file
                        file_info.append(row)
                else:
                    if file.split('.')[-1] == 'txt' and file.split('_')[1][0] == 'B':
                        rna_id = file.split('_')[0]
                        replicate = file.split('_')[1][1:]
                        control_file = [file_name for file_name in files if len(file_name.split('_')) >1 and replicate == file_name.split('_')[1][1:] and file_name != file][0]

                        row = {'rna_id':rna_id, 'date':'', 'experimenter':'','probe':'', 'temperature':'', 'magnesium':'', 'ddNTP':'', 'rna_begin':'',\
                            'rna_end':'','rt_begin_pos':'', 'rt_end_pos':'', 'replicate':replicate, 'probe_file':file, 'control_file':control_file,\
                                'qushape_file':'','reference_qushape_file':'','map_file':'','discard':''}
                        file_info.append(row)
                    elif file.split('.')[-1] in ['qushapey', 'map']:
                        for info in list(data):
                            if file.split('.')[-1] == 'qushapey':
                                row['qushape_file'] = file
                            elif file.split('.')[-1] == 'map':
                                row['map_file'] = file
                            else:
                                row[info] = ''
                        file_info.append(row)

            data = pd.concat([data, pd.DataFrame(file_info)], ignore_index=True)
            data = data.reset_index(drop=True)
            data.index = data.index + 1
            data.index.name = 'id'
            data.to_csv(sample_file, sep='\t', index=True)
        else :
            print("Operation aborted.")


    def refactor(self, action: str = "addpositions"):
        """
        Rename result files to include new conditions.
        """
        extra_config = dict()

        if action == "addpositions":
            extra_config["refactor_addpositions"] = True
            targets = ["all_add_positions"]
        elif action == "rename":
            extra_config["refactor_rename"] = True
            targets = ["all_rename"]

        else:
            raise fire.core.FireError(f"invalid refactor option {action}")

        try:
            sm.snakemake(
                os.path.join(base_path, "workflow", "Snakefile"),
                configfiles=[self._config] if self._config else None,
                config=extra_config,
                targets=targets,
                cores=self._cores,
                keepgoing=self._keepgoing,
                use_conda=True,
                verbose=self._verbose,
                # listrules=True,
                conda_prefix="~/.ipasuite/conda",
            )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)

    def convert_qushape(self, rerun_incomplete=False):
        """
        Convert qushape to qushapey files.
        """
        self._config = self._choose_config(self._config)
        extra_config = dict()
        extra_config["convert_qushape"] = True
        targets = ["all_convert"]

        try:
            sm.snakemake(
                os.path.join(base_path, "workflow", "Snakefile"),
                configfiles=[self._config] if self._config else None,
                config=extra_config,
                targets=targets,
                cores=self._cores,
                keepgoing=self._keepgoing,
                use_conda=True,
                verbose=self._verbose,
                force_incomplete=rerun_incomplete,
                # listrules=True,
                conda_prefix="~/.ipasuite/conda",
            )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)

    def qushape(
        self,
        action="all_reactivity",
        dry_run=False,
    ):
        """
        Open QuShape for each file to treat.
        """
        self._config = self._choose_config(self._config)
        targets = ["all_reactivity"]
        extra_config = dict()

        extra_config["qushape"] = {"run_qushape": True}
        self._cores = 1

        try:
            # if True:
            if self.check():
                sm.snakemake(
                    os.path.join(base_path, "workflow", "Snakefile"),
                    configfiles=[self._config] if self._config else None,
                    config=extra_config,
                    targets=targets,
                    cores=1,
                    keepgoing=self._keepgoing,
                    dryrun=dry_run,
                    use_conda=True,
                    verbose=self._verbose,
                    conda_prefix="~/.ipasuite/conda",
                    rerun_triggers=["mtime"],
                )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)
            logger.error("to get more information, type : ipasuite log")

    def correl(self, dry_run=False, run_qushape=False, rerun_incomplete=False):
        """
        Execute the Ipanemap suite pipeline and stop before the IPANEMAP prediction step.
        """
        self.run(action="all_aggregate",
            dry_run=dry_run,
            run_qushape=run_qushape, 
            rerun_incomplete=rerun_incomplete)
        
    def run(self, action="all", dry_run=False, run_qushape=False, rerun_incomplete=False):
        """
        Execute the Ipanemap suite pipeline.
        """
        self._config = self._choose_config(self._config)
        targets = [action]
        extra_config = dict()

        if run_qushape:
            extra_config["qushape"] = {"run_qushape": True}
            self._cores = 1
        try:
            if self.check():
                # if True:
                sm.snakemake(
                    os.path.join(base_path, "workflow", "Snakefile"),
                    configfiles=[self._config] if self._config else None,
                    config=extra_config,
                    targets=targets,
                    cores=self._cores,
                    keepgoing=self._keepgoing,
                    dryrun=dry_run,
                    use_conda=True,
                    verbose=self._verbose,
                    conda_prefix="~/.ipasuite/conda",
                    rerun_triggers=["mtime"],
                    force_incomplete=rerun_incomplete,
                    #debug_dag = True,
                )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)
            logger.error("to get more information, type : ipasuite log")
#        if os.path.exists('./results/5.2-ipanemap-temp'):
#            shutil.rmtree('./results/5.2-ipanemap-temp')

    def unlock(self):
        """
        Remove Snakemake locks on the working directory.
        
        This runs Snakemake with flag --unlock. It may be necessary to run 
        this after unclean termination. In this case, Snakemake will ask you
        to unlock, but make sure that there is no other still running job.
        """
        self._config = self._choose_config(self._config)
        targets = ["all"]
        extra_config = dict()
        try:
            if self.check():
                # if True:
                sm.snakemake(
                    os.path.join(base_path, "workflow", "Snakefile"),
                    configfiles=[self._config] if self._config else None,
                    config=extra_config,
                    targets=targets,
                    use_conda=True,
                    conda_prefix="~/.ipasuite/conda",
                    rerun_triggers=["mtime"],
                    unlock = True,
                )
        except Exception as e:
            logger.error(traceback.format_exc())
            logger.error(e)
            logger.error("to get more information, type : ipasuite log")


    def log(self, step=None, clean=False, print_filename=False):
        """
        Show pipeline scripts log.
        """
        self._config = self._choose_config(self._config)
        with open(self._config, "r") as cfd:
            config = yaml.load(cfd)
            validate(config, base_path + "/workflow/schemas/config.schema.yaml")
        if step is not None:
            pattern = f"{config['results_dir']}/logs/{step}*.log"
        else:
            pattern = f"{config['results_dir']}/logs/*.log"
        gl = glob(pattern)
        if clean:
            for file in gl:
                os.remove(file)
            if len(gl) > 0:
                logger.info("Log cleaned.")

        else:
            for file in gl:
                with open(file, "r") as fd:
                    read = fd.read()
                    if len(read) > 2:
                        if print_filename:
                            print(file)
                        print(read[:-1])

    def clean(self, from_step="reactivity", keep_log=False):
        """
        Clean pipeline files.
        """
        self._config = self._choose_config(self._config)
        try:
            begin = step_order.index(from_step)
        except ValueError:
            logger.error(f"Authorized value for from_step :{step_order}")
        with open(self._config, "r") as cfd:
            config = yaml.load(cfd)
            validate(config, base_path + "/workflow/schemas/config.schema.yaml")

        for folder in step_order[begin:]:
            try:
                shutil.rmtree(
                    os.path.join(config["results_dir"], config["folders"][folder])
                )
                shutil.rmtree(
                    os.path.join(
                        config["results_dir"], "figures", config["folders"][folder]
                    )
                )
                logger.info(f"{folder} cleaned")
            except FileNotFoundError:
                logger.info(f"no folder for {folder}")

            if not keep_log:
                gl = glob(f"{config['results_dir']}/logs/{folder}*.log")
                for file in gl:
                    os.remove(file)
                if len(gl) > 0:
                    logger.info(f"{folder} logs cleaned")


    def comparison(self):
        """
        Compare secondary structures.
        """
        if 'results' in os.listdir('./'):
            folder = './results/5.2-structure'
            file_names = os.listdir('./results/5.2-structure')
        else:
            folder = './'
            file_names = os.listdir('./')

        def get_sequence_structure(filename):
            with open(f'{filename}', 'r') as file:
                    content = file.read()
            lines = content.split('\n')
            sequence = lines[1]
            structure = lines[2]
            mfe_base_pair = []
            mfe_bp = list(RNA.ptable(structure))[1:]
            for j in range(len(mfe_bp)):
                if mfe_bp[j] != 0:
                    mfe_base_pair.append([j,mfe_bp[j]-1])
            model = {'sequence': sequence, 'secondary_structure':structure, 'mfe_base_pair':mfe_base_pair}
            return model
        
        def model_vs_model(model1, model2):
            identical_bases = []
            different_bps = []
            for n, i in enumerate(model1['secondary_structure']):
                if i == model2['secondary_structure'][n] and model2['secondary_structure'][n] not in ['(',')']:
                    identical_bases.append(1)
                    different_bps.append(1)
                elif i in ['(',')'] and model2['secondary_structure'][n] in ['(',')']:
                    bp = [bp for bp in model1['mfe_base_pair'] if bp[0] == n]
                    if bp[0] in model2['mfe_base_pair']:
                        identical_bases.append(1)
                        different_bps.append(1)
                    else:
                        identical_bases.append(0)
                        different_bps.append(-1)
                elif i == '.' and model2['secondary_structure'][n] == '-':
                    identical_bases.append(1)
                    different_bps.append(1)
                else:
                    identical_bases.append(0)
                    different_bps.append(0)
            model = {'sequence' : model1['sequence'], 'secondary_structure' : model1['secondary_structure'], 'identical_bases':identical_bases, 'different_bps':different_bps}
            return model
        
        def varnaplot(model_info, color, path):
            conda_prefix = os.environ.get("CONDA_PREFIX")
            varnaapi.set_VARNA(f'{conda_prefix}/lib/varna/VARNA.jar')
            v = varnaapi.Structure(sequence=model_info['sequence'], structure=model_info['secondary_structure'])
            v.set_algorithm('radiate')
            v.update(bpStyle="lw", spaceBetweenBases=0.75, bpIncrement=1.3)
            v.add_colormap(values = color, vMin=-1, vMax=1, caption='', style={-1: '#fc0303',-0.9: '#ffffff', 0.9:'#ffffff',1:'#1d05f5'})
            tmpname = path
            v.savefig(tmpname)

        file_list = '\n      '.join(file_names)
        print(f'\nModel DBN filename : \n      {file_list}')
        response1 = input("Please enter DBM filename of the first model : ")
        while response1 not in file_names:
            print(f'\nModel DBN filename : \n      {file_list}')
            response1 = input("Please enter DBM filename of the first model : ")
        
        if response1 in file_names:
            filename1 = f'{folder}/{response1}'
            file_names.remove(response1)
            file_list = '\n      '.join(file_names)
            print(f'\nModel DBN filename : \n      {file_list}')
            response2 = input("Please enter DBM filename of the second model : ")
            while response2 not in file_names:
                print(f'\nModel DBN filename : \n      {file_list}')
                response2 = input("Please enter DBM filename of the second model : ")
            if response2 in file_names:
                filename2 = f'{folder}/{response2}'
                
                model1 = get_sequence_structure(filename1)
                model2 = get_sequence_structure(filename2)
                if model1['sequence'].upper().replace('T','U') != model2['sequence'].upper().replace('T','U') :
                    print('\n*** Please verify the input DBN files, the sequences differ between the two DBN files. ***')

                elif len(model1['secondary_structure']) != len(model2['secondary_structure']):
                    print('\n*** Please verify the input DBN files, the secondary structures in the two DBN files have different lengths. ***')
                
                else:
                    output_name = input("Please enter the output folder name : ")
                    while output_name == '' or output_name == ' ':
                        output_name = input("Please enter the output folder name : ")
                    if output_name != '' and output_name != ' ':
                        output_name = output_name
                        model_vs = model_vs_model(model1, model2)

                        os.makedirs(f'./comparaison_models/{output_name}', exist_ok=True)
                        file_paths = {}
                        for j in model_vs:
                            if j not in ['sequence', 'secondary_structure']:
                                file_paths[f'{j}_file_path'] = f'./comparaison_models/{output_name}/{j}.txt'
                                varnaplot_path = f'./comparaison_models/{output_name}/{j}.varna'
                                varnaplot(model_vs, model_vs[j], varnaplot_path)
            

    def _check_dup(self, path):
        path = pathlib.Path(path)
        if not os.path.exists(path):
            return True
        files = sorted(os.listdir(path))

        distinctFileClass = []

        for file_x in files:
            are_identical = False

            for class_ in distinctFileClass:
                are_identical = filecmp.cmp(
                    path / file_x, path / class_[0], shallow=False
                )
                if are_identical:
                    class_.append(file_x)
                    break
            if not are_identical:
                distinctFileClass.append([file_x])
        contains_dup = False
        for class_ in distinctFileClass:
            if len(class_) > 1:
                contains_dup = True
                logger.warning(
                    "files "
                    + " and ".join(class_)
                    + " are identicals. Make sur your are did not"
                    " make a mistake while preparing you replicates",
                )
        return contains_dup

    def _check_conditions(self, samples, conditions, rna_id, name, type="ipanemap"):
        sample_missing = False
        query = [f' {cname} == "{cond}" &' for cname, cond in conditions.items()] + [
            f' rna_id == "{rna_id}"'
        ]
        query = "".join(query)
        if len(samples.query(query)) == 0:
            sample_missing = True
            logger.error(
                f"{type} pool {name} cannot be handled because "
                f"no sample is available for condition: {query}"
            )
        return sample_missing

    def _check_samples(self, config, samples):
        sample_missing = False
        ipan = config["ipanemap"]["pools"]
        for pool in ipan:
            if "external_conditions" in pool:
                for cond in pool["external_conditions"]:
                    if not os.path.exists(cond["path"]):
                        logger.error(
                            f" ipanemap {cond['path']} not found in "
                            f"{pool['id']} -"
                            f" external {cond['name']}"
                        )
            for conds in pool["conditions"]:
                sample_missing = (
                    self._check_conditions(
                        samples, conds, pool["rna_id"], pool["id"], "ipanemap"
                    )
                    or sample_missing
                )
        for comp in config["footprint"]["compares"]:
            sample_missing = (
                self._check_conditions(
                    samples, comp["condition1"], comp["rna_id"], comp["id"], "footprint"
                )
                or sample_missing
            )
            sample_missing = (
                self._check_conditions(
                    samples, comp["condition2"], comp["rna_id"], comp["id"], "footprint"
                )
                or sample_missing
            )

        return sample_missing

    def check(self, verbose=False):
        """
        Verify input files, sample availability, duplicates, and identical raw files.
        """
        seq_missing = False
        raw_missing = False
        sample_missing = False
        seq_not_fasta = False
        has_duplicate = False
        self._config = self._choose_config(self._config)
        with open(self._config, "r") as cfd:
            config = yaml.load(cfd)
            validate(config, base_path + "/workflow/schemas/config.schema.yaml")

        samples = load_samples.get_unindexed_samples(config)

        indexes = load_samples.get_indexes(config, replicates_in_index=True)
        ind_samples = samples.set_index(indexes, drop=True)
        files = (
            list(samples["probe_file"])
            + list(samples["control_file"])
            + list(samples["qushape_file"])
            + list(samples["reference_qushape_file"])
        )
        files = [
            os.path.join(config["rawdata"]["path_prefix"], str(f))
            for f in files
            if f != "" and not (f != f)
        ]

        seqs = [seq for idx, seq in config["sequences"].items()]

        for file in seqs:
            if not os.path.exists(file):
                seq_missing = True
                logger.error(f"Sequence: {file} not found")
            else:
                try:
                    fastaseqs = fasta.fasta_iter(file)
                    for fname, fseq in fastaseqs:
                        if not fseq.isupper():
                            logger.warning("Sequence contains lowercase")
                        if "T" in fseq:
                            logger.warning(
                                "Sequence contains 'T', it should contains "
                                "'U' in order to be RNA"
                            )
                except Exception as e:
                    logger.error(traceback.format_exc())
                    logger.error(e)
                    logger.error(f"Sequence: {file} does not seems to be a fasta file")
                    seq_not_fasta = True

        for file in files:
            if not os.path.exists(file):
                raw_missing = True
                logger.warning(f"Raw data : {file} not found")

        index_count = pd.Index(ind_samples.index).value_counts()
        for i, val in index_count[index_count > 1].items():
            logger.error(
                f"Found {val} Duplicated samples for"
                f" condition {dict(zip(indexes, i))}"
            )
            logger.error(
                "Every sample of the same condition must has a different replicate id. " 
                f"You must  change replicate id in the samples.tsv file"

            )
            has_duplicate = True

        self._check_dup(
            os.path.join(config["results_dir"], config["folders"]["fluo-fsa"])
        )
        self._check_dup(
            os.path.join(config["results_dir"], config["folders"]["fluo-ceq8000"])
        )

        sample_missing = self._check_samples(config, samples)
        if raw_missing or seq_missing or sample_missing or seq_not_fasta:
            logger.error("-- Problems where found when checking pipeline --")
        else:
            print("Configuration check succeed")

        return not seq_missing and not sample_missing and not seq_not_fasta and not has_duplicate


def main_wrapper():
    fire.Fire(Launcher)
