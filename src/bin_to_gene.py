import logging
from typing import Optional, Union

import numpy as np
from tqdm import tqdm

from .bin import Bin
from .utils import load_gene_file
from .utils import extend_interval

FORMAT = '%(levelname)s: %(message)s'
logging.basicConfig(format=FORMAT)


class BinToGene:
    """
    Given a cell by bin matrix, convert it to a cell by gene matrix.
    This can be useful for e.g. when dealing with sc-ATAC-seq data
    served in cell-by-bin format and one wishes to run analysis on the
    cell-by-gene matrix instead (perhaps due to readily available tools).

    Using a genecode (v34) file, determine the gene location (start, end)
    and extend its extremes by a fixed or gene-length dependent value to
    get a new location interval I. Then sum over all the bins whose coordinates
    intersect with I that belong to the corresponding chromosome.

    For efficiency, binary search is used to determine start and end bins.
    """

    def __init__(
            self,
            gencode_path: str = "gencode.v34.genes.protein_coding.csv",
            operation: str = 'sum',
            extend: Optional[Union[str, int]] = '5x',
            max_extend: Optional[Union[str, int]] = 50000,
            stream_direction: bool = True,
            op_extend: Optional[Union[str, int]] = '1x',
            max_op_extend: Optional[Union[str, int]] = 10000,
            logging_lvl="WARNING"):
        """
        Parameters:
        ___________
        gencode_path: string specifying the path of the gencode file

        operation: Can be 'sum' or 'mean'. Specifies whether to add bins
            that intersect with the gene or take their mean instead.

        extend: String ending in x specifying a multiplier of the gene length
            to add to its extremes, or int specifying exact number of base
            pairs to add.

        max_extend: Same format as extend, but will be used as a threshold.

        stream_direction: Specifies whether to take stream direction
            into consideration. If set to True, then op_stream_extend specifies
            the value that will be added to the opposite stream end.

        op_extend: Same format as extend, but applied to the opposite
            stream end only. Ignored if stream_direction set to False.

        max_op_extend: Same as max_extend applied to op_stream_extend.

        logging_lvl: Specify logging level.
        """
        self.logger = logging.getLogger('BinToGene')
        self.logger.setLevel(getattr(logging, logging_lvl))

        self.gencode = load_gene_file(gencode_path, self.logger)
        self.extend = extend
        self.max_extend = max_extend
        self.stream_direction = stream_direction
        self.op_extend = op_extend
        self.max_op_extend = max_op_extend

        if operation == 'sum':
            self.op = np.sum
        else:
            self.op = np.mean

    def binary_search(
            self,
            bin_list: np.ndarray,
            start: int,
            end: int,
            search_first: bool = True):
        """
        Run binary search to determine first or last bin that intersects with
        the given (start, end) range.

        Parameters:
        ___________
        bin_list: Array of Bin.

        start, end: Ints specifying the interval to search for.

        search_first: If True, look for the first bin that intersects
            the given range. If False, look for the last bin instead.

        Returns:
        ________
        result: Integer specifying the index of the first or last bin
            in bin_list. Will return -1 if no such bin is found.
        """
        i = 0
        j = len(bin_list) - 1
        result = -1

        while i <= j:
            m = (i + j) // 2
            if bin_list[m].suceeds(start, end):
                j = m - 1
            elif bin_list[m].precedes(start, end):
                i = m + 1
            else:
                result = m
                if search_first:
                    j = m - 1
                else:
                    i = m + 1

        return result

    def get_gene_counts(self, gene, bin_dict, x):
        """
        Given a gene, find the intersecting bins in bin_dict and add
        up the corresponding columns of x to return a vector of counts
        for the given gene.

        Parameters:
        ___________
        gene: Single row DataFrame

        bin_dict: Dictionary of bins. Keys are chromosome names.

        x: Cell by bin matrix

        Returns:
        ________
        counts: 1D np.ndarray, or None if no intersection found

        """
        start, end = gene['start'], gene['end']
        seqname = gene['seqname']

        if seqname not in bin_dict:
            return None

        # Determine whether different extension should be used for
        # different stream direction
        if self.stream_direction is None or self.stream_direction == False:
            strand = None
        else:
            strand = gene['strand']

        # Extend the gene interval accordingly
        start, end = extend_interval(
            start, end, self.extend, self.max_extend,
            strand, self.op_extend, self.max_op_extend)

        # Perform binary search to determine the first and last bins
        # of the same sequence that intersect with the gene range

        # WARNING: These are the first and last indices in the
        # appropriate list bin_dict[seqname] and not in the
        # original data matrix
        first = self.binary_search(bin_dict[seqname], start, end, True)
        last = self.binary_search(bin_dict[seqname], start, end, False)

        assert first <= last, "First bin greater index than last."

        # Add up all the intersecting bins
        if first != -1 and last != -1:  # intersection found
            true_start_bin = bin_dict[seqname][first].index
            true_end_bin = bin_dict[seqname][last].index
            gene_counts = self.op(x[:, true_start_bin:true_end_bin + 1], axis=1)

            return gene_counts

        return None  # if no intersection found

    def convert(
            self,
            x: np.ndarray,
            bin_names: Union[np.ndarray, list],
            prefix: Optional[str] = '',
            delim1: str = ':',
            delim2: str = '-',
            to_dense: bool = False):
        """
        Given a cell by bin matrix convert it to cell by gene.

        x: Cell by bin matrix. Can be a numpy array or sparse matrix.

        bin_names: Names of bins/columns.

        prefix: If the names contain a prefix such as in 'chr1'. The remanining
            chromosome names should be 1,...,22 and X or Y.

        delim1: Delimiter that separates sequence name from range.

        delim2: Delimiter that separates start from end.
        """
        if x.shape[1] != len(bin_names):
            raise ValueError("Mismatch between number of bin and length"
                             "of bin names.")

        if prefix is None:
            prefix = ""

        # Dictionary of chromosomes
        bin_dict = {'chr' + str(i): [] for i in list(range(1, 23)) + ['X', 'Y']}

        for i, bin_name in enumerate(bin_names):
            # Split bin name into sequence name and interval
            try:
                seqname, interval = str(bin_name).split(delim1)
                start, end = interval.split(delim2)
                start, end = int(start), int(end)
            except:
                raise ValueError("Invalid bin name encountered.")

            seqname = 'chr' + seqname[len(prefix):]  # Remove prefix and add chr
            if seqname in bin_dict:  # Only accept chr1-chr22 and chrX, chrY
                bin_dict[seqname].append(Bin(seqname, start, end, i))

        # Will be using binary search, so we sort bins by start value
        for seq in bin_dict:
            bin_dict[seq].sort(key=lambda x: x.start)
            bin_dict[seq] = np.array(bin_dict[seq])

        active_gene_ids = []
        counts = []

        for index in tqdm(self.gencode.index):
            gene = self.gencode.loc[index]
            gene_counts = self.get_gene_counts(gene, bin_dict, x)
            if gene_counts is not None:
                counts.append(gene_counts)
                active_gene_ids.append(gene['gene_id'])

        return np.array(counts), np.array(active_gene_ids)


if __name__ == '__main__':
    import anndata
    a = anndata.read_h5ad('src/cell_by_bin.h5ad')
    b = BinToGene()
    counts, names = b.convert(a.X, a.var_names)
