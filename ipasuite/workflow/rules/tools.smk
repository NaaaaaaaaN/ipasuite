

if RAW_DATA_TYPE == "fluo-ceq8000":
    rule fluo_ceq8000:
        input: construct_path(step="fluo-ceq8000")
        output: protected(construct_path(step="fluo-ce"))
        log: construct_path('fluo-ce', ext=".log", log_dir=True)
        message: f"Converting ceq8000 data for qushape: {MESSAGE} replicate"
                 f" {{wildcards.replicate}}"
        shell:
            f"ceq8000_to_tsv {{input}} {{output}} &> {{log}}"

if config["qushape"]["use_subsequence"]:
    rule split_fasta:
        input:
            ancient(get_refseq)
        output:
            f"{RESULTS_DIR}/{config['folders']['subseq']}/{{rna_id}}_{{rt_end_pos}}-{{rt_begin_pos}}.fasta"
    
        message:
            f"Fragmenting fasta file {{input}} "
            f"from {{wildcards.rt_end_pos}} to {{wildcards.rt_begin_pos}}"
    
        log:
            f"logs/{config['folders']['subseq']}/{{rna_id}}_{{rt_end_pos}}-{{rt_begin_pos}}.log"
        shell:
            f"split_fasta {{input}} {{output}} --begin "
            f"{{wildcards.rt_end_pos}}"
            f" --end {{wildcards.rt_begin_pos}}"

if "fluo-ce" in RAW_DATA_TYPE:
    rule generate_project_qushape:
        input:
            rx = ancient(construct_path("fluo-ce")),
            bg = ancient(construct_path("fluo-ce", control = True )),
            refseq = ancient(lambda wildcards: get_subseq(wildcards, split_seq=True)),
            refproj = ancient(get_qushape_refproj)
        message: f"Generate QuShape project for {MESSAGE}"
                 f"- replicate {{wildcards.replicate}}"
        params:
            refseq=lambda wildcards, input: f"--refseq={input.refseq}",
            refproj=lambda wildcards, input: expand('--refproj={refproj}', refproj=input.refproj)[0] if len(input.refproj) > 0 else "",
            ddNTP=lambda wildcards: f"--ddNTP={get_ddntp_qushape(wildcards)}",
            channels= construct_dict_param(config["qushape"], "channels"),
            overwrite="--overwrite=untreated"
        log: construct_path('qushape', ext=".log", log_dir=True, split_seq=True)
        output: construct_path("qushape", ext=".qushapey", split_seq=True)
        shell:
            f"qushape_proj_generator {{input.rx}} {{input.bg}}"
            f" {{params}} --output={{output}} &> {{log}}"

if "fluo-fsa" in RAW_DATA_TYPE:
    rule generate_project_qushape:
        input:
            rx = ancient(construct_path("fluo-fsa", ext=".fsa")),
            bg = ancient(construct_path("fluo-fsa", control = True, ext=".fsa" )),
            refseq = ancient(lambda wildcards: get_subseq(wildcards, split_seq=True)),
            refproj = ancient(get_qushape_refproj)
        message: f"Generate QuShape project for {MESSAGE}"
                 f"- replicate {{wildcards.replicate}}"
        params:
            refseq=lambda wildcards, input: f"--refseq={input.refseq}",
            refproj=lambda wildcards, input: expand('--refproj={refproj}', refproj=input.refproj)[0] if len(input.refproj) > 0 else "",
            ddNTP=lambda wildcards: f"--ddNTP={get_ddntp_qushape(wildcards)}",
            channels= construct_dict_param(config["qushape"], "channels"),
            overwrite="--overwrite=untreated"
        log: construct_path('qushape', ext=".log", log_dir=True, split_seq=True)
        output: construct_path("qushape", ext=".qushapey", split_seq=True)
        shell:
            f"qushape_proj_generator {{input.rx}} {{input.bg}}"
            f" {{params}} --output={{output}} &> {{log}}"



rule extract_reactivity:
    """
    Extract reactivities from qushapey files
    
    Writes reactivity.tsv report, which directly contains the data from qushape's report.
    Moreover, it produces a simple bar plot of the data (areas RX BG, and their difference) for each sequence position
    """
    input:
        qushape = construct_path("qushape", ext=".qushapey", split_seq=True),
        refseq = ancient(lambda wildcards: get_subseq(wildcards, split_seq=True))
    output:
        react=construct_path("reactivity", split_seq=True),
        plot=report(construct_path("reactivity", ext=".reactivity.svg",  figure=True, split_seq=True), category="3.1-Reactivity", subcategory=CONDITION),
        #,protect = protected(construct_path("qushape", ext=".qushape"))
    params:
        rna_file=lambda wildcards, input: f"--rna_file {input.refseq}" if config["qushape"]["check_integrity"] else "",
        plot_title=lambda wildcards: f"--plot_title='Reactivity of {MESSAGE.format(wildcards=wildcards)} - replicate {wildcards.replicate}'",
        run_qushape=lambda wildcards: f"--launch-qushape --qushape-conda-env={config['qushape']['run_qushape_conda_env']}" if config["qushape"]["run_qushape"] else "",
    message: f"Extracting reactivity from QuShape for {MESSAGE}"
             f"- replicate {{wildcards.replicate}}"
    log: construct_path('reactivity', ext=".log", log_dir=True, split_seq=True)

    shell:
        f"qushape_extract_reactivity {{input.qushape}} {{params}}"
        f" --output={{output.react}} --plot={{output.plot}} &> {{log}}"

rule normalize_reactivity:
    """
    Normalizes the 'raw' reactivities, given as areas, in the reactivity.tsv files

    See rule extract_reactivity: these areas come directly from qushape

    The script normalize reactivity produces a tsv file normreact.tsv (main output) as well as map and shape files
    Moreover, it produces profile plots of position-wise normalized reactivities 
    (comparing two ways of normalization: /simple/ and /interquartile/.)  
    
    Note: the shape and map files start at positions 1, even if the first reactivities are not contained in files reactivity.tsv.
    For this purpose, the first positions of undefined reactivities are padded with 'undefined' entries in the normalize_reactivity script.
    """
    input: construct_path("reactivity", split_seq=True)
    output:
        nreact=construct_path("normreact", split_seq=True),
        plot=report(construct_path("normreact", ext=".normreact.svg",
            figure=True, split_seq=True) ,
                category="3.2-Normalized reactivity", subcategory=CONDITION) if
        config["normalization"]["plot"] else [],
        shape_file = construct_path("normreact", ext=".shape", split_seq=True),
        map_file = construct_path("normreact", ext=".map", split_seq=True),
    message: f"Normalizing reactivity for {MESSAGE}"
             f" - replicate {{wildcards.replicate}}"
    log: construct_path('normreact', ext=".log", log_dir=True, split_seq=True)
    params:
        reactive_nucleotides = get_reactive_nucleotides,
        stop_percentile = construct_param(CNORM, "stop_percentile"),
        low_norm_reactivity_threshold = construct_param(CNORM, "low_norm_reactivity_threshold"),
        norm_methods = construct_list_param(CNORM, "norm_methods"),
        shape_file_normcol = construct_normcol(),
        snorm_out_perc= construct_param(CNORM, "simple_outlier_percentile"),
        simple_norm_term_avg_percentile = construct_param(CNORM, "simple_norm_term_avg_percentile"),
        plot = lambda wildcards, output: f"--plot={output.plot}" if config["normalization"]["plot"] else "",
        plot_title=lambda wildcards: f"--plot_title='Normalized Reactivity of {MESSAGE.format(wildcards=wildcards)} - replicate {wildcards.replicate}'"
    shell:
        f"normalize_reactivity {{params}} {{input}}"
        f" --shape_output={{output.shape_file}}"
        f" --map_output={{output.map_file}}"
        f" --output={{output.nreact}}  &> {{log}}"

if config["qushape"]["use_subsequence"]:
    rule align_reactivity_to_ref:
        input: unpack(get_align_reactivity_inputs)
        output: construct_path("alignnormreact")
        params:
            rt_end_pos = get_align_begin,
            #rna_end = get_align_end

        log: construct_path('alignnormreact', ext=".log", log_dir=True)
        shell:
            f"shift_reactivity {{input.norm}} {{input.refseq}} {{output}}"
            f" --begin {{params.rt_end_pos}} &> {{log}}"
            #f" --end {{params.rna_end}} &> {{log}}"


rule aggregate_reactivity:
    input:
        norm= lambda wildcards: expand(construct_path(aggregate_input_type()),
                replicate=get_replicate_list(wildcards), allow_missing=True),
        refseq = get_refseq
    output:
        full= construct_path("aggreact", show_replicate = False),
        shape_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".shape"),
        shape_IP_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".ip.shape"),
        map_file = construct_path("aggreact-ipanemap", show_replicate=False, ext=".map"),
        relation_file = construct_path("aggreact-ipanemap", show_replicate=False, ext="_consistency.csv"),
        plot =report(construct_path("aggreact", ext=".aggreact.svg",
            show_replicate=False, figure=True),
            category="4-Aggregated reactivity", subcategory=CONDITION) if config["aggregate"]["plot"] else [],
        fullplot = report(construct_path("aggreact", ext=".aggreact.full.svg",
            show_replicate=False, figure=True),
            category="4-Aggregated reactivity", subcategory=CONDITION) if config["aggregate"]["plot"] else [],

    #message: f"Aggregating normalized reactivity for {MESSAGE}"
    log: construct_path('aggreact', ext=".log", log_dir=True, show_replicate=False)
    params:
        norm_method= construct_normcol(),
        min_std = construct_param(config["aggregate"], "min_std"),
        reactivity_medium = construct_param(config["aggregate"], "reactivity_medium"),
        reactivity_high = construct_param(config["aggregate"], "reactivity_high"),
#        minndp = construct_param(config["aggregate"], "min_ndata_perc"),
#        mindndp = construct_param(config["aggregate"], "min_nsubdata_perc"),
#        maxmp = construct_param(config["aggregate"], "max_mean_perc"),
#        mind = construct_param(config["aggregate"], "min_dispersion"),
        plot = lambda wildcards, output: f"--plot={output.plot} --fullplot={output.fullplot}" if
        config["aggregate"]["plot"] else "",
        plot_title=lambda wildcards: f"--plot_title='Average reactivity of {MESSAGE.format(wildcards=wildcards)}'",
        refseq = lambda wildcards, input: (expand('--refseq={refseq}', refseq=input.refseq)[0] if len(input.refseq) > 0 else ""),
    shell:
        f"aggregate_reactivity {{input.norm}}"
        f" --output={{output.full}} {{params}}"
        f" --shape_output={{output.shape_file}}"
        f" --shape_IP_output={{output.shape_IP_file}}"
        f" --map_output={{output.map_file}}"
        f" --relation_output={{output.relation_file}}"
        f" --err_on_dup={config['aggregate']['err_on_dup']} &> {{log}}"

rule footprint:
    input:
        # the two aggreact.tsv files that are compared by this footprint
        tsv = get_footprint_inputs,
        
        # structure as dbn to be annotated by footprinting information 
        struct_1 = lambda wildcards: expand(
            f"{RESULTS_DIR}/{{folder}}/{wildcards.rna_id}_pool_{get_footprint_pool_id(wildcards)}_1.dbn",
            folder=config["folders"]["structure"],
            allow_missing=True,
        ),
        struct_2 = lambda wildcards: expand(
            f"{RESULTS_DIR}/{{folder}}/{wildcards.rna_id}_pool_{get_footprint_pool_id(wildcards)}_2.dbn",
            folder=config["folders"]["structure"],
            allow_missing=True,
        ),

    output:
        tsv=f"{RESULTS_DIR}/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}.tsv",
        plot=report(f"{RESULTS_DIR}/figures/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}.svg", category="5-Footprint", subcategory="{rna_id} - {foot_id}"),
        diff_plot=report(f"{RESULTS_DIR}/figures/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}_difference.svg", category="5-Footprint", subcategory="{rna_id} - {foot_id}"),
        structure_plot_1_svg=report(f"{RESULTS_DIR}/figures/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}_structure_1.svg", category="5-Footprint", subcategory="{rna_id} - {foot_id}"),
        structure_plot_1_varna=report(f"{RESULTS_DIR}/figures/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}_structure_1.varna", category="5-Footprint", subcategory="{rna_id} - {foot_id}"),
        structure_plot_2_svg=report(f"{RESULTS_DIR}/figures/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}_structure_2.svg", category="5-Footprint", subcategory="{rna_id} - {foot_id}"),
        structure_plot_2_varna=report(f"{RESULTS_DIR}/figures/{config['folders']['footprint']}/{{rna_id}}_footprint_{{foot_id}}_structure_2.varna", category="5-Footprint", subcategory="{rna_id} - {foot_id}"),
    log: "results/logs/footprint_{rna_id}_footprint_{foot_id}.log"
    params:
        ttest_pvalue_thres = construct_param(config["footprint"]["config"],
                "ttest_pvalue_thres"),
   #     deviation_type = construct_param(config["footprint"]["config"],
   #             "deviation_type"),
        diff_thres = construct_param(config["footprint"]["config"],
                "diff_thres"),
        ratio_thres = construct_param(config["footprint"]["config"],
                "ratio_thres"),
        cond1_name = lambda wildcards: f"--cond1_name='{get_footprint_condition_name(wildcards, 1)}'",
        cond2_name = lambda wildcards: f"--cond2_name='{get_footprint_condition_name(wildcards, 2)}'",
        plot_title = lambda wildcards: f"--plot_title='Compared reactivity between {get_footprint_condition_name(wildcards, 1)} and {get_footprint_condition_name(wildcards, 2)}'",
        diff_plot_title = lambda wildcards: f"--diff_plot_title='Difference between {get_footprint_condition_name(wildcards, 1)} and {get_footprint_condition_name(wildcards, 2)}'"
    shell:
        "footprint {input.tsv} --output={output.tsv} {params}"
        " --plot={output.plot} --diff_plot={output.diff_plot} --plot_format=svg; "
        "footprint {input.tsv} {params} --structure={input.struct_1} --structure_plot={output.structure_plot_1_svg}; "
        "footprint {input.tsv} {params} --structure={input.struct_1} --structure_plot={output.structure_plot_1_varna}; "
        "footprint {input.tsv} {params} --structure={input.struct_2} --structure_plot={output.structure_plot_2_svg}; "
        "footprint {input.tsv} {params} --structure={input.struct_2} --structure_plot={output.structure_plot_2_varna}"
