import numpy as np
import cloudpickle

__deletion = 'deletion'
__DaMP = 'damp'

class GIData(object):
    def __init__(self, 
                 values, 
                 rows, 
                 cols,
                 row_meta=None,
                 col_meta=None,
                 has_duplicates=False,
                 check_symmetric=True):
        assert all([gene.isupper() for gene in rows])
        assert all([gene.isupper() for gene in cols])
        
        assert isinstance(values, np.ndarray), type(values)
        assert values.dtype == float, values.dtype

        assert isinstance(rows, np.ndarray)
        assert rows.dtype.type == np.str_, rows.dtype.type

        assert isinstance(cols, np.ndarray)
        assert cols.dtype.type == np.str_, cols.dtype.type

        assert len(rows.shape) == 1
        assert len(cols.shape) == 1

        self.values = values
        self.rows = rows
        self.cols = cols
        self.n_rows = len(rows)
        self.n_cols = len(cols)
        self.shape = self.values.shape

        assert(self.values.shape == (self.n_rows, self.n_cols))

        if has_duplicates:
            assert self.no_duplicates(self.rows)
            assert self.no_duplicates(self.cols)
        
        self.has_duplicates = has_duplicates
        self.is_symmetric = check_symmetric and np.allclose(self.values, self.values.T, equal_nan=True)

    
    @staticmethod
    def no_duplicates(genes):
        return len(genes) == len(set(genes))
    
    @classmethod
    def from_multiIndexDF(cls, mi_df, gene_name_key = 'Gene'):
        # TODO: add DF metadata
        gis = cls(values = mi_df.values,
                  rows = mi_df.index.get_level_values(gene_name_key).values.astype(str),
                  cols = mi_df.index.get_level_values(gene_name_key).values.astype(str))
        return gis
    
    def to_dict(self):
        return dict(
            values = self.values,
            rows   = self.rows,
            cols   = self.cols,
            has_duplicates = self.has_duplicates
        )
    
    @classmethod
    def from_dict(cls, d):
        return cls(**d)

    def save(self, fp):
        with open (fp, 'wb') as f:
            cloudpickle.dump(self.to_dict(), f)
    
    @classmethod
    def load(cls, fp):
        with open (fp, 'rb') as f:
            return cls.from_dict(cloudpickle.load(f))
