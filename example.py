from src import BinToGene
import anndata

a = anndata.read_h5ad('src/resources/cell_by_bin.h5ad')[:100]
btg = BinToGene(n_jobs=4)
counts, ids = btg.convert(a.X, a.var_names.to_numpy())

