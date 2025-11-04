"""
Quality control and filtering rules
"""

rule filter_cells:
    """
    Filter low-quality cells based on QC metrics
    - min/max genes detected
    - mitochondrial content
    - doublet detection
    """
    input:
        h5ad = lambda wildcards: datasets_config[wildcards.dataset]["input_path"]
    output:
        h5ad = "{results}/{dataset}/filtered.h5ad",
        qc_plots = directory("{reports}/{dataset}/qc_plots")
    params:
        min_genes = lambda wildcards: datasets_config[wildcards.dataset].get("min_genes", 200),
        max_genes = lambda wildcards: datasets_config[wildcards.dataset].get("max_genes", 6000),
        max_pct_mt = lambda wildcards: datasets_config[wildcards.dataset].get("max_pct_mt", 15),
        min_cells_per_gene = lambda wildcards: datasets_config[wildcards.dataset].get("min_cells_per_gene", 3)
    log:
        "{logs}/qc/{dataset}_filter.log"
    conda:
        "../../environment.yml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 8000 * attempt,
        runtime = 120
    script:
        "../../scripts/filter_normalize.py"

rule qc_report:
    """
    Generate comprehensive QC HTML report
    """
    input:
        raw = lambda wildcards: datasets_config[wildcards.dataset]["input_path"],
        filtered = "{results}/{dataset}/filtered.h5ad"
    output:
        report = "{reports}/{dataset}/qc_report.html"
    log:
        "{logs}/qc/{dataset}_report.log"
    conda:
        "../../environment.yml"
    threads: 2
    resources:
        mem_mb = 4000
    script:
        "../../scripts/generate_qc_report.py"