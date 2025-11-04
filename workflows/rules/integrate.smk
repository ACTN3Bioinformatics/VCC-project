"""
Batch integration rules (optional)
"""

rule batch_correction:
    """
    Correct batch effects using Harmony or BBKNN
    Only runs if batch_correction: true in config
    """
    input:
        "results/{dataset}/balanced.h5ad"
    output:
        "results/{dataset}/integrated.h5ad"
    params:
        batch_key = lambda wildcards: config['datasets'][wildcards.dataset].get("batch_key", "batch"),
        method = lambda wildcards: config['datasets'][wildcards.dataset].get("integration_method", "harmony"),
        run = lambda wildcards: config['datasets'][wildcards.dataset].get("batch_correction", False)
    log:
        "logs/integrate/{dataset}.log"
    conda:
        str(Path(workflow.basedir) / "environment.yml")
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 120
    script:
        str(Path(workflow.basedir) / "scripts" / "integration.py")