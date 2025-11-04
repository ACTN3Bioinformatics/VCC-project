"""
Feature engineering rules
Extract biological features (pathways, TFs, etc.)
"""

rule extract_features:
    """
    Extract biological features:
    - Gene pathway membership
    - TF targets
    - Regulatory networks
    - H1-hESC specific features
    """
    input:
        h5ad = lambda wildcards: (
            "{results}/{dataset}/integrated.h5ad" 
            if datasets_config[wildcards.dataset].get("batch_correction", False)
            else "{results}/{dataset}/balanced.h5ad"
        )
    output:
        h5ad = "{results}/{dataset}/features.h5ad",
        feature_list = "{results}/{dataset}/feature_list.txt"
    params:
        use_pathways = True,
        use_tf_targets = True,
        use_regulatory_networks = True,
        databases = ["KEGG", "Reactome", "GO_Biological_Process"]
    log:
        "{logs}/features/{dataset}.log"
    conda:
        "../../environment.yml"
    threads: 4
    resources:
        mem_mb = 12000,
        runtime = 90
    script:
        "../../scripts/feature_engineering.py"

rule combine_features_splits:
    """
    Combine features with train/val/test splits
    """
    input:
        features = "{results}/{dataset}/features.h5ad",
        train = "{results}/{dataset}/splits/train.h5ad",
        val = "{results}/{dataset}/splits/val.h5ad",
        test = "{results}/{dataset}/splits/test.h5ad"
    output:
        final = "{results}/{dataset}/final.h5ad"
    log:
        "{logs}/features/{dataset}_combine.log"
    conda:
        "../../environment.yml"
    threads: 2
    script:
        "../../scripts/combine_features.py"