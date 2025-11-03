"""
Data download rules
Downloads and prepares demonstration datasets from Zenodo/Figshare

Source: Replogle et al. 2022 K562 Essential Perturb-seq
- Primary: Zenodo (scPerturb database)
- Backup: Figshare (original authors) 
- No login required, direct download
"""

rule download_demo_data:
    """
    Download Replogle et al. 2022 K562 essential dataset
    
    Data source: 
    - Zenodo record 7041849 (scPerturb curated, PRIMARY)
    - Figshare (original authors, BACKUP)
    
    Original: 1.44 GB, ~188k cells, ~2000 perturbations
    Output: ~500MB, ~10k cells, ~150 perturbations (subsetted)
    
    Processing time: ~5-10 minutes (download + subset)
    Memory required: ~4GB RAM during processing
    """
    output:
        h5ad = "data_local/demo/replogle_subset.h5ad",
        metadata = "data_local/demo/metadata.json"
    params:
        # PRIMARY: Zenodo (scPerturb database) - TESTED AND WORKING
        url = config.get("demo_data_url", 
                        "https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_essential.h5ad"),
        # BACKUP: Figshare (original authors) - TESTED AND WORKING
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
        runtime = 30  # 30 minutes timeout
    script:
        "../../scripts/download_demo_data.py"

rule verify_data:
    """
    Verify integrity of downloaded data
    - Check file can be loaded
    - Validate structure (cells, genes, perturbations)
    - Confirm QC metrics exist
    """
    input:
        "data_local/demo/replogle_subset.h5ad"
    output:
        touch("data_local/demo/.verified")
    log:
        "logs/download/verify_demo.log"
    conda:
        "../../environment.yml"
    script:
        "../../scripts/verify_data.py"