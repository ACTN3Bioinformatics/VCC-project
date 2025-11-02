"""
Feature engineering rules
"""

rule extract_features:
    input:
        h5ad = lambda wildcards: (
            f"{{results}}/{{dataset}}/integrated.h5ad" 
            if datasets_config[wildcards.dataset].get("batch_correction", False)
            else f"{{results}}/{{dataset}}/balanced.h5ad"
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
        "logs/features/{dataset}.log"
    conda:
        "../../environment.yml"
    threads: 4
    resources:
        mem_mb = 12000,
        runtime = 90
    script:
        "../../scripts/feature_engineering.py"

rule combine_features_splits:
    input:
        features = "{results}/{dataset}/features.h5ad",
        train = "{results}/{dataset}/splits/train.h5ad",
        val = "{results}/{dataset}/splits/val.h5ad",
        test = "{results}/{dataset}/splits/test.h5ad"
    output:
        final = "{results}/{dataset}/final.h5ad"
    log:
        "logs/features/{dataset}_combine.log"
    conda:
        "../../environment.yml"
    threads: 2
    script:
        "../../scripts/combine_features.py"