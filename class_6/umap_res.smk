# perform clustering
rule umap_res:
    input:
        get_data_path
    output:
        umap_res_results = os.path.join(result_path,'{analysis}','res_0.6.rds'),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","umap_res_results_{analysis}.log"),
    params:
        partition=config.get("partition"),
        assay = lambda w: annot_dict["{}".format(w.analysis)]["assay"],
        metadata = lambda w: annot_dict["{}".format(w.analysis)]["metadata"],
        control = lambda w: annot_dict["{}".format(w.analysis)]["control"],
        logfc_threshold = config["logfc_threshold"],
        test_use = config["test_use"],
        min_pct = config["min_pct"],
        return_thresh = config["return_thresh"],
        #resolution = config["resolution"],
    script:
        "../scripts/umap_res.R"