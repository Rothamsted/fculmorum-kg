from rpy2.robjects.vectors import StrVector
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
import pandas as pd
#import numpy as np
import sys

def r2pd_dataframe(df):
    """ Convert R dataframe to Pandas DataFrame
        Requires an R DataFrame
        Outputs a Pandas DataFrame """
    with localconverter(ro.default_converter + pandas2ri.converter):
      pd_from_r_df = ro.conversion.rpy2py(df)
    return pd_from_r_df

def get_fungi_mart(ds):
    return R.useMart(biomart = 'fungi_mart',
                     host = "https://fungi.ensembl.org/",
                     dataset=ds)

def return_bm_df(atts, mart):
    atts_vector = StrVector([x for x in atts])
    BM = R.getBM(attributes = atts_vector, mart = mart)
    return r2pd_dataframe(BM)

base = sys.argv[1]

R.library("biomaRt")

print(R.listMarts(host = "https://fungi.ensembl.org/"))

# Get the fungi marts
mart_names = ['fculmorum_eg_gene' ,'fgraminearum_eg_gene', 'fpseudograminearum_eg_gene']
marts = [get_fungi_mart(mart) for mart in mart_names]
f_culmorum_mart, f_graminearum_mart, f_pseudogram_mart = marts[0], marts[1], marts[2]

# Get the attributes and the filters
f_culmorum_attributes, f_graminearum_attributes, f_pseudogram_attriibtes = R.listAttributes(f_culmorum_mart), R.listAttributes(f_graminearum_mart), R.listAttributes(f_pseudogram_mart)

f_cul_att_df = r2pd_dataframe(f_culmorum_attributes)
gram_df = f_cul_att_df.name[f_cul_att_df.name.str.startswith('f')]
gram_df.to_string()

f_cul_att_df[f_cul_att_df['name'].str.contains("gene")]

f_culmorum_atts = ['ensembl_gene_id', 'external_gene_name', 'fpseudograminearum_eg_homolog_ensembl_gene',
                    'fpseudograminearum_eg_homolog_associated_gene_name']
f_cul_pd = return_bm_df(f_culmorum_atts, f_culmorum_mart)
# Remove empty columns
f_cul_pd_update = f_cul_pd[f_cul_pd.fpseudograminearum_eg_homolog_ensembl_gene != '']

f_coding_type_atts = ['ensembl_gene_id', 'external_gene_name', 'gene_biotype']
f_cul_coding_pd = return_bm_df(f_coding_type_atts, f_culmorum_mart)
f_cul_coding_pd = f_cul_coding_pd[f_cul_coding_pd.gene_biotype == 'protein_coding']

merged_df = f_cul_coding_pd.merge(f_cul_pd_update, how = 'inner', on=['ensembl_gene_id', 'ensembl_gene_id'])
merged_df = merged_df[merged_df.gene_biotype == 'protein_coding']
del merged_df['external_gene_name_x']
merged_df = merged_df[['ensembl_gene_id', 'external_gene_name_y', 'gene_biotype', 'fpseudograminearum_eg_homolog_ensembl_gene', 'fpseudograminearum_eg_homolog_associated_gene_name']]
merged_df

merged_df.to_csv(f'{base}/homologs_f_cul.txt', sep="\t", index=None, header=True)
