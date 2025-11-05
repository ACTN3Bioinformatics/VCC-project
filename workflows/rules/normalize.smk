"""
Normalization and scaling rules
"""

rule normalize:
    """
    Normalize counts and log-transform
    - Total count normalization
    - Log1p transformation
    - Optional: regress out technical variables
    """
    input:
        "results/{dataset}/filtered.h5ad"
    output:
        "results/{dataset}/normalized.h5ad"
    params:
        target_sum = lambda wildcards: config['datasets'][wildcards.dataset].get("target_sum", 10000),
        log_transform = lambda wildcards: config['datasets'][wildcards.dataset].get("log_transform", True),
        regress_out = lambda wildcards: config['datasets'][wildcards.dataset].get("regress_out", [])
    log:
        "logs/normalize/{dataset}.log"
    conda:
        str(Path(workflow.basedir) / "environment.yml")
    threads: 4
    resources:
        mem_mb = 12000,
        runtime = 60
    script:
        str(Path(workflow.basedir) / "scripts" / "filter_normalize.py")


rule scale:
    """
    Scale genes to unit variance and zero mean
    """
    input:
        "results/{dataset}/normalized.h5ad"
    output:
        "results/{dataset}/scaled.h5ad"
    params:
        scale = lambda wildcards: config['datasets'][wildcards.dataset].get("scale", True),
        max_value = 10
    log:
        "logs/scale/{dataset}.log"
    conda:
        str(Path(workflow.basedir) / "environment.yml")
    threads: 4
    resources:
        mem_mb = 12000,
        runtime = 30
    script:
        str(Path(workflow.basedir) / "scripts" / "filter_normalize.py")