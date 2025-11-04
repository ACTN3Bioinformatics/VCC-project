"""
Batch integration rules (optional)
"""

rule batch_correction:
    """
    Correct batch effects using Harmony or BBKNN
    Only runs if batch_correction: true in config
    """
    input:
        "{results}/{dataset}/balanced.h5ad"
    output:
        "{results}/{dataset}/integrated.h5ad"
    params:
        batch_key = lambda wildcards: datasets_config[wildcards.dataset].get("batch_key", "batch"),
        method = lambda wildcards: datasets_config[wildcards.dataset].get("integration_method", "harmony"),
        run = lambda wildcards: datasets_config[wildcards.dataset].get("batch_correction", False)
    log:
        "{logs}/integrate/{dataset}.log"
    conda:
        "../../environment.yml"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 120
    script:
        "../../scripts/integration.py"