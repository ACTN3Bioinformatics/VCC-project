"""
Data download rules
"""

rule download_demo_data:
    output:
        h5ad = "{data_local}/demo/replogle_subset.h5ad",
        metadata = "{data_local}/demo/metadata.json"
    params:
        url = config.get("demo_data_url", 
                        "https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_essential.h5ad"),
        alternative_url = config.get("demo_alternative_url",
                                    "https://plus.figshare.com/ndownloader/files/42444534"),
        max_cells = config.get("demo_max_cells", 10000),
        target_perturbations = config.get("demo_target_perturbations", 150)
    log:
        "logs/download/demo_data.log"
    conda:
        "../../environment.yml"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 30
    script:
        "../../scripts/download_demo_data.py"

rule verify_data:
    input:
        "{data_local}/demo/replogle_subset.h5ad"
    output:
        touch("{data_local}/demo/.verified")
    log:
        "logs/download/verify_demo.log"
    conda:
        "../../environment.yml"
    script:
        "../../scripts/verify_data.py"