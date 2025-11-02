"""
Unit tests for QC and filtering functions
"""

import pytest
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


@pytest.fixture
def mock_adata():
    """Create mock AnnData object for testing"""
    np.random.seed(42)
    n_obs = 1000
    n_vars = 500
    
    X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars))
    obs = pd.DataFrame({
        'n_genes_by_counts': np.random.randint(100, 5000, n_obs),
        'total_counts': np.random.randint(1000, 50000, n_obs),
        'pct_counts_mt': np.random.uniform(0, 30, n_obs),
        'target_gene': np.random.choice(['GENE_A', 'GENE_B', 'GENE_C'], n_obs)
    })
    var = pd.DataFrame(index=[f'Gene_{i}' for i in range(n_vars)])
    
    return AnnData(X=X, obs=obs, var=var)


def test_filter_cells_min_genes(mock_adata):
    """Test filtering cells by minimum genes"""
    min_genes = 200
    adata_filtered = mock_adata[mock_adata.obs['n_genes_by_counts'] >= min_genes]
    assert (adata_filtered.obs['n_genes_by_counts'] >= min_genes).all()
    assert adata_filtered.n_obs <= mock_adata.n_obs


def test_filter_cells_max_mt(mock_adata):
    """Test filtering cells by max mitochondrial content"""
    max_pct_mt = 15
    adata_filtered = mock_adata[mock_adata.obs['pct_counts_mt'] <= max_pct_mt]
    assert (adata_filtered.obs['pct_counts_mt'] <= max_pct_mt).all()
    assert adata_filtered.n_obs <= mock_adata.n_obs


def test_normalization(mock_adata):
    """Test count normalization"""
    target_sum = 10000
    sc.pp.normalize_total(mock_adata, target_sum=target_sum)
    sums = np.array(mock_adata.X.sum(axis=1)).flatten()
    assert np.allclose(sums, target_sum, rtol=0.01)


def test_log_transform(mock_adata):
    """Test log transformation"""
    sc.pp.log1p(mock_adata)
    assert (mock_adata.X >= 0).all()
    assert mock_adata.X.max() < 20


if __name__ == "__main__":
    pytest.main([__file__, "-v"])