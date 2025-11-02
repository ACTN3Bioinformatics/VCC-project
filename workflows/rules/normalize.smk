"""
Normalization rules
"""

rule normalize:
    input:
        "{results}/{dataset}/filtered.h5ad"
    output:
        "{results}/{dataset}/normalized.h5ad"
    params:
        target_sum = lambda wildcards: datasets_config[wildcards.dataset].get("target_sum", 10000),
        log_transform = lambda wildcards: datasets_config[wildcards.dataset].get("log_transform", True),
        regress_out = lambda wildcards: datasets_config[wildcards.dataset].get("regress_out", ["total_counts", "pct_counts_mt"])
    log:
        "logs/normalize/{dataset}.log"
    conda:
        "../../environment.yml"
    threads: 4
    resources:
        mem_mb = 12000,
        runtime = 60
    script:
        "../../scripts/filter_normalize.py"

rule scale:
    input:
        "{results}/{dataset}/normalized.h5ad"
    output:
        "{results}/{dataset}/scaled.h5ad"
    params:
        scale = lambda wildcards: datasets_config[wildcards.dataset].get("scale", True),
        max_value = 10
    log:
        "logs/scale/{dataset}.log"
    conda:
        "../../environment.yml"
    threads: 4
    resources:
        mem_mb = 12000,
        runtime = 30
    script:
        "../../scripts/filter_normalize.py"