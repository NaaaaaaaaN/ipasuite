TOOLS = "workflow/scripts/tools/"


if RAW_DATA_TYPE == "fluo-ceq8000":
    rule fluo_ceq8000:
        conda: "../envs/tools.yml"
        input: construct_path(step="fluo-ceq8000", results_dir=False)
        output: protected(construct_path(step="fluo-ce", results_dir=False))
        log: construct_path('fluo-ce', ext=".log", log_dir=True)
        message: f"Converting ceq8000 data for qushape: {MESSAGE} replicate"
                 f" {{wildcards.replicate}}"
        shell:
            f"python {TOOLS}/ceq8000_to_tsv.py {{input}} {{output}} &> {{log}}"



# If a qushape file from outside exists
#rule import_qushape:
#    input: construct_path("qushape", ext=".qushape", results_dir=False)
#    output: construct_path("qushape", ext=".qushape")
#    log: construct_path('qushape', ext=".log", log_dir=True)
#    message: f"Importing from ressource: {MESSAGE} replicate"
#             f"{{wildcards.replicate}}"
#    shell:
#        "cp {input} {output} &> {log}"

if config["qushape"]["use_subsequence"]:
    rule split_fasta:
        conda: "../envs/tools.yml"
        input:
            ancient(get_refseq)
        output:
            f"results/{config['folders']['subseq']}/{{rna_id}}_{{rt_end_pos}}-{{rt_begin_pos}}.fasta"
    
        message:
            f"Fragmenting fasta file {{input}} "
            f"from {{wildcards.rt_end_pos}} to {{wildcards.rt_begin_pos}}"
    
        log:
            "logs/{config['folders']['subseq']}/{rna_id}_{rt_end_pos}-{rt_begin_pos}.log"
        shell:
            f"python {TOOLS}/split_fasta.py {{input}} {{output}} --begin "
            f"{{wildcards.rt_end_pos}}"
            f" --end {{wildcards.rt_begin_pos}}"

rule generate_project_qushape:
    conda: "../envs/tools.yml"
    input:
        rx = ancient(construct_path("fluo-ce", results_dir = False)),
        bg = ancient(construct_path("fluo-ce", control = True, results_dir = False)),
        refseq = ancient(lambda wildcards: get_subseq(wildcards, split_seq=True)),
        refproj = ancient(get_qushape_refproj)
    message: f"Generate QuShape project for {MESSAGE}"
             f"- replicate {{wildcards.replicate}}"
    params:
        refseq=lambda wildcards, input: expand('--refseq={refseq}', refseq=input.refseq)[0] if len(input.refseq) > 0 else "",
        refproj=lambda wildcards, input: expand('--refproj={refproj}', refproj=input.refproj)[0] if len(input.refproj) > 0 else "",
        ddNTP= construct_param(config["qushape"], "ddNTP"),
        channels= construct_list_param(config["qushape"], "channels")
        #TODO channels
        #channels=
    #output: protected(construct_path("qushape", ext=".qushape"))

    log: construct_path('qushape', ext=".log", log_dir=True, split_seq=True)
    output: construct_path("qushape", ext=".qushape", split_seq=True)
    shell:
        f"python {TOOLS}/qushape_proj_generator.py {{input.rx}} {{input.bg}}"
        f" {{params}} --output={{output}} &> {{log}}"

rule extract_reactivity:
    conda:  "../envs/tools.yml"
    input: construct_path("qushape", ext=".qushape", split_seq=True)
    output:
        react=construct_path("reactivity", split_seq=True),
        plot=report(construct_path("reactivity", ext=".reactivity.svg", 
            figure=True, split_seq=True),
            category="3.1-Reactivity", subcategory=CONDITION)
        #,protect = protected(construct_path("qushape", ext=".qushape"))
    message: f"Extracting reactivity from QuShape for {MESSAGE}"
             f"- replicate {{wildcards.replicate}}"
    log: construct_path('reactivity', ext=".log", log_dir=True, split_seq=True)

    shell:
        f"""
        set +e
        python {TOOLS}/qushape_extract_reactivity.py {{input}} \
                --output={{output.react}} --plot={{output.plot}} &> {{log}}

        exitcode=$?
        if [ $exitcode != 0 ] ; then
            echo -n 'ERROR: '
            cat {{log}}
        fi
        exit $exitcode
        """

rule normalize_reactivity:
    conda:  "../envs/tools.yml"
    input: construct_path("reactivity", split_seq=True)
    output:
        nreact=construct_path("normreact", split_seq=True),
        plot=report(construct_path("normreact", ext=".normreact.svg",
            figure=True, split_seq=True),
                category="3.2-Normalized reactivity", subcategory=CONDITION)
    message: f"Normalizing reactivity for {MESSAGE}"
             f" - replicate {{wildcards.replicate}}"
    log: construct_path('normreact', ext=".log", log_dir=True, split_seq=True)
    params:
        react_nuc = construct_list_param(CNORM, "reactive_nucleotides"),
        st_perc = construct_param(CNORM, "stop_percentile"),
        low_norm_reac_thres = construct_param(CNORM, "low_norm_reactivity_threshold"),
        norm_methods = construct_list_param(CNORM, "norm_methods"),
        snorm_out_perc= construct_param(CNORM, "simple_outlier_percentile"),
        snorm_term_avg_perc= construct_param(CNORM, "simple_norm_term_avg_percentile")
    shell:
        f"python {TOOLS}/normalize_reactivity.py {{params}} {{input}}"
        f" --output={{output.nreact}} --plot={{output.plot}} &> {{log}}"

if config["qushape"]["use_subsequence"]:
    rule align_reactivity_to_ref:
        conda: "../envs/tools.yml"
        input: unpack(get_align_reactivity_inputs)
        output: construct_path("alignnormreact")
        params:
            rt_end_pos = get_align_begin,
            #rna_end = get_align_end

        log: construct_path('alignnormreact', ext=".log", log_dir=True)
        shell:
            f"python {TOOLS}/shift_reactivity.py {{input.norm}} {{input.refseq}}"
            f" {{output}} --begin {{params.rt_end_pos}} &> {{log}}"
            #f" --end {{params.rna_end}} &> {{log}}"
     


rule aggregate_reactivity:
    conda:  "../envs/tools.yml"
    input:
        norm= lambda wildcards: expand(construct_path(aggregate_input_type()),
                replicate=get_replicate_list(wildcards), allow_missing=True),
        refseq = lambda wildcards: get_refseq(wildcards)
    output:
        full= construct_path("aggreact", show_replicate = False),
        compact = construct_path("aggreact-ipanemap", show_replicate=False,
                ext=".txt"),
        plot =report(construct_path("aggreact", ext=".aggreact.svg",
            show_replicate=False, figure=True),
            category="4-Aggregated reactivity", subcategory=CONDITION),
        fullplot =report(construct_path("aggreact", ext=".aggreact.full.svg",
            show_replicate=False, figure=True),
            category="4-Aggregated reactivity", subcategory=CONDITION)

    #message: f"Aggregating normalized reactivity for {MESSAGE}"
    log: construct_path('aggreact', ext=".log", log_dir=True, show_replicate=False)
    params:
        norm_method= construct_normcol(),
        minndp = construct_param(config["aggregate"], "min_ndata_perc"),
        mindndp = construct_param(config["aggregate"], "min_nsubdata_perc"),
        maxmp = construct_param(config["aggregate"], "max_mean_perc"),
        mind = construct_param(config["aggregate"], "min_dispersion"),
        #refseq = lambda wildcards, input: expand('--refseq={refseq}', refseq=input.refseq)[0] if len(input.refseq) > 0 else ""
    shell:
        f"python {TOOLS}/aggregate_reactivity.py {{input.norm}}"
        f" --output={{output.full}} {{params}}"
        f" --ipanemap_output={{output.compact}}"
        f" --plot={{output.plot}} --fullplot={{output.fullplot}} &> {{log}}"

#rule ipanemap:
#    conda: "../envs/ipanemap.yml"
#    input: construct_path("aggreact-ipanemap", replicate = False)
#    output: "python workflow/scripts/IPANEMAP/IPANEMAP.py"
#
