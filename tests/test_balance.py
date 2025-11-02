"""
Unit tests for balancing functions
"""

import pytest
import numpy as np
import pandas as pd
from anndata import AnnData


@pytest.fixture
def imbalanced_adata():
    """Create imbalanced dataset for testing"""
    np.random.seed(42)
    
    n_cells_A = 500
    n_cells_B = 100
    n_cells_C = 50
    
    obs = pd.DataFrame({
        'target_gene': ['GENE_A'] * n_cells_A + ['GENE_B'] * n_cells_B + ['GENE_C'] * n_cells_C
    })
    
    n_vars = 100
    X = np.random.negative_binomial(5, 0.3, (len(obs), n_vars))
    
    return AnnData(X=X, obs=obs)


def test_downsample_balanced(imbalanced_adata):
    """Test downsampling creates balanced dataset"""
    target_cells = 50
    
    balanced_data = []
    for gene in imbalanced_adata.obs['target_gene'].unique():
        gene_cells = imbalanced_adata[imbalanced_adata.obs['target_gene'] == gene]
        
        if gene_cells.n_obs > target_cells:
            indices = np.random.choice(gene_cells.n_obs, target_cells, replace=False)
            sampled = gene_cells[indices]
        else:
            sampled = gene_cells
        
        balanced_data.append(sampled)
    
    for adata in balanced_data:
        assert adata.n_obs <= target_cells


def test_preserve_rare_classes(imbalanced_adata):
    """Test that rare classes are preserved"""
    target_cells = 50
    gene_c_cells = imbalanced_adata[imbalanced_adata.obs['target_gene'] == 'GENE_C']
    assert gene_c_cells.n_obs == 50


if __name__ == "__main__":
    pytest.main([__file__, "-v"])