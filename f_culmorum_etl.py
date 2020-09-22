import numpy as np
import requests
import wget
import gzip, shutil
import os
import pandas as pd
from Bio import SeqIO
import urllib
import json
import sys
from optparse import OptionParser
from sh import gunzip

def mkfile(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def gunzip(x):
    print("Unzipping now...\n")
    with gzip.open(x, 'rb') as f_in:
        with open(x.replace('.gz', ''), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print("Unzipped!")

def unzip_tidy(fn, path):
    gz_f = lambda x: gunzip(x) if '.gz' in x else None
    gz_f(fn)
    return fn.replace('.gz', '')

def blast_extract(path, fn, rm):
    """ Cuts the first 2 columns (gene and proteins) and returns the dataframe """
    df = pd.read_csv(f'{path}/{fn}', sep="\t", header=None )
    df = df[[0,1]]
    df.columns = ['Gene', 'Protein']
    df = df.drop_duplicates()
    if rm:
        os.remove(os.path.join(f'{path}', fn))
    return df

def blast_concat(paths, fn):
    """
    Path is a list of the separate paths, fns is the static filename
    Concats all the BLAST outputs into a bigger mapping file
    """
    df_list = []
    for path in paths:
        print(f"Parsing {path}{fn}")
        try:
            df = blast_extract(path, fn, False)
            df_list.append(df)
        except:
            print("failed")
    return df_list

def df_split_col_delimiter(df, target_column, separator):
    """ df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split

    Returns: a dataframe with each entry for the target column separated, with each element moved into a new row.
    The values in the other columns are duplicated across the newly divided rows.
    """
    def splitListToRows(row, row_accumulator, target_column, separator):
        split_row = row[target_column].split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(splitListToRows, axis=1, args=(
        new_rows, target_column, separator))
    new_df = pd.DataFrame(new_rows)
    return new_df

def concat_dropna(list):
    mapping_df = pd.concat(list)
    return mapping_df.dropna()

def request_data(url, fn_out):
    r = requests.get(url, stream=True, headers={'User-agent': 'Mozilla/5.0'})
    if r.status_code == 200:
        with open(fn_out, 'wb+') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
    return fn_out

def split_df(df, col, delim, index):
    return df[col].str.split(delim).str[index]

def mapping_pd(df, df_type , df_tax_id, name, list, o):

        species_mapping = { 'query': df['query_name'],
                                name: df_type,
                                'tax_id': df_tax_id
                                }

        mapping_df = pd.DataFrame(species_mapping, columns=['query', name, 'tax_id' ])
        mapping_df['species_name'] = o
        list.append(mapping_df)

def ko_desc_mapping(df, df_2, df_tax_id, list, o):

        species_mapping = { 'query': df['query_name'],
                                'KO': df_2['KEGG_ko'],
                                'description': df_2['eggNOG free text description'],
                                'tax_id': df_tax_id
                                }

        mapping_df = pd.DataFrame(species_mapping, columns=['query', 'KO', 'description', 'tax_id' ])
        mapping_df['species_name'] = o
        list.append(mapping_df)
        return mapping_df

def taxid_rest(o):

    ebi_url = f"https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/any-name/{o}"
    r = urllib.request.urlopen(ebi_url)
    data = json.loads(r.read().decode(r.info().get_param('charset') or 'utf-8'))
    return data[0]['taxId']

def ensembl_data_url(gen_ver, species, dir):

    # Get species specific url
    gene_url, peptide_url = ensembl_url(species, gen_ver)
    # Make file if it doesn't exist
    mkfile(f'{dir}/{species}/ensembl')
    # DL datasets
    gene_fn = wget.download(gene_url, f'{dir}/{species}/ensembl')
    peptide_fn = wget.download(peptide_url, f'{dir}/{species}/ensembl')
    # unzip
    gz_f = lambda x: gunzip(x) if '.gz' in x else None
    gz_f(gene_fn), gz_f(peptide_fn)
    files = os.listdir(f'{dir}/{species}/ensembl/')
    #remove gunzip files
    [os.remove(os.path.join(f'{dir}/{species}/ensembl/', file)) for file in files if file.endswith(".gz")]
    return gene_fn.replace('.gz', ''), peptide_fn.replace('.gz', '')


def ensembl_url(species, gen_ver):
    # These can be extended with more species
    gene_url = {
                'fusarium_culmorum': f'ftp://ftp.ensemblgenomes.org/pub/release-{gen_ver}/fungi/gff3/fusarium_culmorum/Fusarium_culmorum.EF1.{gen_ver}.gff3.gz'
                }[species]

    peptide_url = {
                    'fusarium_culmorum': f'ftp://ftp.ensemblgenomes.org/pub/fungi/release-{gen_ver}/fasta/fusarium_culmorum/pep/Fusarium_culmorum.EF1.pep.all.fa.gz'
                    }[species]

    return gene_url, peptide_url

def yield_records(dir, char):
    """
    Yields all relevant fasta records and replaces the ID to be uniprot only. Can be modified as needed.
    Returns the updated generator, avoids pointless iterations.
    """
    with open(dir) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            seq_record.id = seq_record.id.split(char)[1]
            seq_record.id = seq_record.description = seq_record.id.replace('.seq','')
            yield seq_record

def yield_ensembl_records(x):
    with open(x) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            seq_record.description = ""
            seq_record.id = seq_record.id.replace("P", "T")
            yield seq_record

def ensembl_records(dir, char):
    """
    Yields all relevant fasta records and replaces the ID to be uniprot only. Can be modified as needed.
    Returns the updated generator, avoids pointless iterations.
    """
    with open(dir) as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            seq_record.id = seq_record.id.split(char)[1]
            seq_record.id = seq_record.description = seq_record.id.replace('.seq','')
            yield seq_record


def uniprot_ftp(species_dict, species):

    species_name = list(species_dict.keys())
    if "fusarium" in species_name:
        # Download
        f_culmorum_uniprot_fn = wget.download('https://www.uniprot.org/uniprot/?query=taxonomy:5516&format=xml&force=true&compress=yes', f'{base}/uniprot/f_culmorum.xml.gz')
        f_culmorum_uniprot_names_fn = wget.download('http://www.uniprot.org/uniprot/?query=taxonomy:5516+AND+reviewed:yes&columns=id,genes(PREFERRED)&format=tab&compress=no', f'{base}/uniprot/fcul_gene_names.txt') # Fusarium Cul proteins
        fungi_all_names_fn = wget.download('http://www.uniprot.org/uniprot/?query=taxonomy:4890+AND+reviewed:yes&columns=id,genes(PREFERRED)&format=tab&compress=no', f'{base}/uniprot/fungi_all_gene_names.txt') # Ascomycota Proteins
        all_uniprot_xml = wget.download('https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:4890&format=xml&force=true&compress=yes', f'{base}/uniprot/all_uniprot.xml.gz') # Ascomycota Proteins
        all_fasta_proteins = wget.download('https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:4890&format=fasta&force=true&compress=no', f'{base}/uniprot/all_fungi_proteins.fa')
        # Unzip
        unzip_tidy(f_culmorum_uniprot_fn, f'{base}/uniprot'), unzip_tidy(f_culmorum_uniprot_names_fn, f'{base}/uniprot'),  unzip_tidy(all_uniprot_xml, f'{base}/uniprot')
        # Rename
        os.rename(f_culmorum_uniprot_fn, f"{base}/uniprot/fculmorum_uniprot.xml.gz"), os.rename(all_uniprot_xml, f"{base}/uniprot/all_fungal_proteins_uniprot.xml.gz")

    for i, s in enumerate(species_name):
        print(f"Writing out to {species_name[i]}")
        fn = wget.download(f'https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:{species_dict[s]}&format=fasta&force=true&compress=no', f'{species_name[i]}/{s}_proteins.fa')
        # Tidy the headers up
        with open(f'{species_name[i]}/{s}_uniprot_proteins.fa', "w+") as handle:
            SeqIO.write(yield_records(f'{species_name[i]}/{s}_proteins.fa', '|'), handle, "fasta")
        os.remove(f'{species_name[i]}/{s}_proteins.fa')
    print("Finished all UniProt downloads")

def retrieve_fasta_files():

    # Fusarium oxysporum
    fn_f_oxy_ensembl = wget.download('ftp://ftp.ensemblgenomes.org/pub/fungi/release-47/fasta/fusarium_oxysporum/pep/Fusarium_oxysporum.FO2.pep.all.fa.gz', f'{base}/ensembl/fusarium_oxysporum.pep.fa.gz')
    fn_f_oxy_ensembl = unzip_tidy(fn_f_oxy_ensembl, f'{base}/ensembl/')
    # Read fasta file then filter it
    with open(f'{fn_f_oxy_ensembl}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_f_oxy_ensembl), handle, "fasta")

    # A nidulans agdb
    fn_f_nidulans_agdb = wget.download('http://www.aspergillusgenome.org/download/sequence/A_nidulans_FGSC_A4/current/A_nidulans_FGSC_A4_current_orf_coding.fasta.gz', f'{base}/agdb/a_nidulans_orf.fa.gz')
    fn_f_nidulans_agdb = unzip_tidy(fn_f_nidulans_agdb, '')

    with open(f'{fn_f_nidulans_agdb}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_f_nidulans_agdb), handle, "fasta")

    # A nidulans Ensembl
    fn_a_nid_ensembl = wget.download('ftp://ftp.ensemblgenomes.org/pub/fungi/release-47/fasta/aspergillus_nidulans/pep/Aspergillus_nidulans.ASM1142v1.pep.all.fa.gz', f'{base}/ensembl/a_nidulans.pep.fa.gz')
    fn_a_nid_ensembl = unzip_tidy(fn_a_nid_ensembl, '')
    with open(f'{fn_a_nid_ensembl}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_a_nid_ensembl), handle, "fasta")

    # A nidulans UniProt
    fn_a_nidulans_uniprot = wget.download('https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:162425&format=fasta&force=true', f'{base}/uniprot/a_nidulans_uniprot.fa')
    with open(fn_a_nidulans_uniprot, "r+") as handle:
        SeqIO.write(yield_records(fn_a_nidulans_uniprot, '|'), handle, "fasta")

    # Fusarium graminearum
    fn_gram_ensembl = wget.download('ftp://ftp.ensemblgenomes.org/pub/fungi/release-47/fasta/fusarium_graminearum/pep/Fusarium_graminearum.RR1.pep.all.fa.gz', f'{base}/ensembl/fusarium_graminearum.pep.fa.gz')
    fn_gram_ensembl = unzip_tidy(fn_gram_ensembl, f'{base}/ensembl/')
    with open(f'{fn_gram_ensembl}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_gram_ensembl), handle, "fasta")

    # Fusarium gram UniProt
    fn_gram_uniprot = wget.download('https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:5506&format=fasta&force=true', f'{base}/uniprot/fusarium_gram_uniprot.fa')
    with open(fn_gram_uniprot, "r+") as handle:
        SeqIO.write(yield_records(fn_gram_uniprot, '|'), handle, "fasta")

    # Fusarium lang Ensembl
    fn_lang_ensembl = wget.download('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/fusarium_langsethiae_gca_001292635/pep/Fusarium_langsethiae_gca_001292635.ASM129263v1.pep.all.fa.gz', f'{base}/ensembl/fusarium_lang_uniprot.fa.gz')
    fn_lang_ensembl = unzip_tidy(fn_lang_ensembl, f'')
    with open(f'{fn_lang_ensembl}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_lang_ensembl), handle, "fasta")

    # Fusarium pesudo Ensembl
    fn_psuedo_ensembl = wget.download('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fusarium_pseudograminearum/pep/Fusarium_pseudograminearum.GCA_000303195.1.pep.all.fa.gz', f'{base}/ensembl/fusarium_pseudogram.pep.fa.gz')
    fn_psuedo_ensembl = unzip_tidy(fn_psuedo_ensembl, '')
    with open(f'{fn_psuedo_ensembl}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_psuedo_ensembl), handle, "fasta")

    # Fusarium pesudo uniprot
    fn_pesudo_uniprot = wget.download('https://www.uniprot.org/uniprot/?query=proteome:UP000007978&format=fasta&force=true', f'{base}/uniprot/fusarium_pseudogram_uniprot.fa')
    with open(fn_pesudo_uniprot, "r+") as handle:
        SeqIO.write(yield_records(fn_pesudo_uniprot, '|'), handle, "fasta")

    # Fusarium venenatum Ensembl
    fn_venea_ensembl = wget.download('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota4_collection/fusarium_venenatum_gca_900007375/pep/Fusarium_venenatum_gca_900007375.ASM90000737v1.pep.all.fa.gz', f'{base}/ensembl/fusarim_venenatum.fa.gz')
    fn_venea_ensembl = unzip_tidy(fn_venea_ensembl, '')
    with open(f'{fn_venea_ensembl}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_venea_ensembl), handle, "fasta")

    # Fusarium evenen UniProt
    fn_venea_uniprot = wget.download('https://www.uniprot.org/uniprot/?query=taxonomy:56646&format=fasta&force=true', f'{base}/uniprot/fusarium_evenen_uniprot.fa')
    with open(fn_venea_uniprot, "r+") as handle:
        SeqIO.write(yield_records(fn_venea_uniprot, '|'), handle, "fasta")

    # Magna Oryzae Ensembl
    fn_magna_oryzae_ensembl = wget.download('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/magnaporthe_oryzae/pep/Magnaporthe_oryzae.MG8.pep.all.fa.gz', f'{base}/ensembl/magna_oryzae_ensembl.pep.all.fa.gz')
    fn_magna_oryzae_ensembl = unzip_tidy(fn_magna_oryzae_ensembl, '')
    with open(f'{fn_magna_oryzae_ensembl}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_magna_oryzae_ensembl), handle, "fasta")

    # Magna Oryzae UniProt
    fn_magna_oryzae_uniprot = wget.download('https://www.uniprot.org/uniprot/?query=magnaporthe%20oryzae&format=fasta&force=true&sort=score&fil=reviewed:yes', f'{base}/uniprot/magna_oryzae_uniprot.fa')
    with open(fn_magna_oryzae_uniprot, "r+") as handle:
        SeqIO.write(yield_records(fn_magna_oryzae_uniprot, '|'), handle, "fasta")

    # Ncrassa UniProt
    fn_ncrassa_uniprot = wget.download('https://www.uniprot.org/uniprot/?query=neurospora%20crassa&format=fasta&force=true&sort=score&fil=reviewed:yes', f'{base}/uniprot/ncrassa_uniprot.fa')
    with open(fn_ncrassa_uniprot, "r+") as handle:
        SeqIO.write(yield_records(fn_ncrassa_uniprot, '|'), handle, "fasta")

    # Secrev UniProt
    fn_s_cerevisiae_uniprot = wget.download('https://www.uniprot.org/uniprot/?query=saccharomyces%20cerevisiae&format=fasta&force=true&sort=score&fil=reviewed:yes', f'{base}/uniprot/s_cerevisiae_uniprot.fa')
    with open(fn_s_cerevisiae_uniprot, "r+") as handle:
        SeqIO.write(yield_records(fn_s_cerevisiae_uniprot, '|'), handle, "fasta")

    # Secrev YGD
    fn_s_cerevisiae_YGD = wget.download('http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz', f'{base}/YGD/s_cerevisiae_YGD.fa.gz')
    fn_s_cerevisiae_YGD = unzip_tidy(fn_s_cerevisiae_YGD, '')
    ygd_id, hgnc = [], []
    with open(fn_s_cerevisiae_YGD, "r+") as handle:
        for v in handle:
            if ">" in v:
                ygd_id.append(v.split(" ")[0].replace(">", ""))
                hgnc.append(v.split(" ")[1])
    ygd_mapping_df = pd.DataFrame({'YGD ID': ygd_id, 'HGNC': hgnc})
    ygd_mapping_df.to_csv(f'{base}/mapping/ygd_hgnc_mapping.txt', sep="\t", index=None)
    with open(f'{fn_s_cerevisiae_YGD}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_s_cerevisiae_YGD), handle, "fasta")

    # Zymo Ensembl
    fn_z_trici = wget.download('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/zymoseptoria_tritici/pep/Zymoseptoria_tritici.MG2.pep.all.fa.gz', f'{base}/ensembl/zymoseptoria_tritici.fa.gz')
    fn_z_trici = unzip_tidy(fn_z_trici, '')
    with open(f'{fn_z_trici}', "r+") as handle:
        SeqIO.write(yield_ensembl_records(fn_z_trici), handle, "fasta")

    # Zymo UniProt
    fn_z_trici_uniprot = wget.download('https://www.uniprot.org/uniprot/?query=zymoseptoria&format=fasta&force=true&sort=score&fil=organism:%22Zymoseptoria%20tritici%20ST99CH_1A5%20[1276529]%22', f'{base}/uniprot/zymoseptoria_tritici_uniprot.fa')
    with open(fn_z_trici_uniprot, "r+") as handle:
        SeqIO.write(yield_records(fn_z_trici_uniprot, '|'), handle, "fasta")

def egg_nog(path, fn_egg, fn_fasta):
    """
        Edits the eggNog paper and returns all the relavant files necessary
    """
    main_df = pd.read_csv(f'{path}/{fn_egg}', sep="\t", error_bad_lines=False, header=None)
    main_df.columns = ['query_name', 'seed eggNOG ortholog', 'seed ortholog evalue', 'seed ortholog score',
                        'Predicted taxonomic group', 'Predicted protein name', 'Gene Ontology terms', 'EC number',
                        'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC',
                        'CAZy', 'BiGG Reaction', 'tax_scope', 'eggNOG OGs', 'bestOG', 'COG Functional Category',
                        'eggNOG free text description']

    othors = main_df['Predicted taxonomic group'].unique()
    othor_pd_dict, othor_taxid = {}, {}
    mapping_pro_list, mapping_gene_list, KEGG_list, EC_list, desc_list = [], [], [], [], []
    # Import the fasta file as a dictionary with the key being the protein name
    fasta_dict = SeqIO.to_dict(SeqIO.parse(open(f'{path}/{fn_fasta}'), "fasta"))
    mapping_path = path.replace("eggNog", "mapping")
    fasta_path = path.replace("eggNog", "BLAST")
    # Fetch based on species
    for o in othors:
        othor_name = o.replace(" ", "_").lower()
        df = main_df.loc[main_df['Predicted taxonomic group'] == o]
        df_tax_id = taxid_rest(othor_name.split("_")[0])
        othor_taxid[othor_name] = df_tax_id
        kegg_desc = df[['query_name', 'KEGG_ko', 'eggNOG free text description']].dropna()

        mapping_pd(df=df, df_type = df['seed eggNOG ortholog'].str.split(".").str[-1], df_tax_id=df_tax_id, name = 'protein_match', list=mapping_pro_list, o=o)
        mapping_pd(df=df, df_type = df['Predicted protein name'], df_tax_id=df_tax_id, name = 'gene_name', list=mapping_gene_list, o=o) # N.B. this is the gene name, not protein.
        mapping_pd(df=df, df_type=df['EC number'], df_tax_id=df_tax_id, name="EC", list=EC_list, o=o)
        kegg_desc['species_name'] = o
        KEGG_list.append(kegg_desc)

        mkfile(f'{fasta_path}/species/{othor_name}')
        with open(f"{fasta_path}/species/{othor_name}/fcul_{othor_name}.fa", "w") as handle:
            for x in df['query_name']:
                SeqIO.write(fasta_dict[x], handle, "fasta")
        print(f'Finished processing eggNog data for species {o}\n')

    # Protein
    mapping_pro_df = concat_dropna(mapping_pro_list)
    mapping_pro_df = mapping_pro_df[~mapping_pro_df['protein_match'].apply(lambda x: len(x) <= 1)] # Remove these non gene names!
    mapping_pro_df['query'] = [x.replace("T", "G") for x in mapping_pro_df['query']]
    mapping_pro_df.to_csv(f'{mapping_path}/fcul_pro_mapping.txt', sep="\t", index=None)
    # Gene
    mapping_gene_df = concat_dropna(mapping_gene_list)
    mapping_gene_ids = [x.replace("T", "G") for x in mapping_gene_df['query']]
    mapping_gene_df['query'] = mapping_gene_ids
    #print(mapping_gene_ids)
    mapping_gene_df['gene_name'] = mapping_gene_df['gene_name'].str.upper()
    mapping_gene_df = mapping_gene_df.drop_duplicates('gene_name', keep='first')
    mapping_gene_df.to_csv(f'{mapping_path}/fcul_gene_mapping.txt', sep="\t", index=None)
    print("Finished writing mapping file")

    # KEGG
    kegg_ortho_df = concat_dropna(KEGG_list)
    kegg_ortho_df = df_split_col_delimiter(df=kegg_ortho_df, target_column='KEGG_ko', separator=',')
    kegg_ortho_df = kegg_ortho_df[kegg_ortho_df['KEGG_ko'].str.contains(':')]
    kegg_ortho_df.to_csv(f'{mapping_path}/fcul_kegg.txt', sep="\t", index=None)

    # EC
    EC_df = concat_dropna(EC_list)
    print(EC_df.columns)
    EC_df = EC_df[EC_df['EC'].str.contains('.')]
    EC_df.to_csv(f'{mapping_path}/fcul_EC.txt', sep="\t", index=None)
    print(f"Written file out to {mapping_path}/fcul_EC.txt")

    return othor_taxid

def string_ppi_data(base):

    #Foxy
    fn_foxy_ppi_string = request_data(url='https://stringdb-static.org/download/protein.links.v11.0/5507.protein.links.v11.0.txt.gz', fn_out=f"{base}/string/Fusarium_oxysporum_ppi.txt.gz")
    fn_foxy_ppi_string = unzip_tidy(fn_foxy_ppi_string,'')
    foxy_ppi_df = pd.read_csv(fn_foxy_ppi_string, sep=" ")
    foxy_ppi_df['protein1'], foxy_ppi_df['protein2'] = split_df(foxy_ppi_df, 'protein1', '.', 1), split_df(foxy_ppi_df, 'protein2', '.', 1)
    foxy_ppi_df['protein1'], foxy_ppi_df['protein2'] = split_df(foxy_ppi_df, 'protein1', 'P0', 0), split_df(foxy_ppi_df, 'protein2', 'P0', 0)
    foxy_ppi_df.to_csv(f'{base}/string/foxy_ppi_stringdb.txt', sep="\t", index=None)
    # foxy atts
    fn_foxy_attributes = request_data(url='https://stringdb-static.org/download/protein.info.v11.0/5507.protein.info.v11.0.txt.gz', fn_out=f'{base}/string/fusarium_oxy_atts.txt.gz')
    fn_foxy_attributes = unzip_tidy(fn_foxy_attributes, '')
    foxy_attributes_df = pd.read_csv(fn_foxy_attributes, sep="\t")
    foxy_attributes_df['protein_external_id']  = split_df(foxy_attributes_df, 'protein_external_id', '.', 1)
    foxy_attributes_df['protein_external_id'] = split_df(foxy_attributes_df, 'protein_external_id', 'P0', 0)
    foxy_attributes_df.to_csv(f'{base}/string/fusarium_oxy_attributes.txt', sep="\t", index=None)

    # NEED TO MERGE and repeat for all others
    # F gram
    fn_fgram_ppi_string = request_data(url='https://stringdb-static.org/download/protein.links.v11.0/5518.protein.links.v11.0.txt.gz', fn_out=f"{base}/string/Fusarium_gram_ppi.txt.gz")
    fn_fgram_ppi_string = unzip_tidy(fn_fgram_ppi_string, '')
    fgram_ppi_df = pd.read_csv(fn_fgram_ppi_string, sep=" ")
    fgram_ppi_df['protein1'], fgram_ppi_df['protein2'] = split_df(fgram_ppi_df, 'protein1', '.', 1), split_df(fgram_ppi_df, 'protein2', '.', 1)

    fgram_ppi_df.to_csv(f'{base}/string/fgram_ppi_stringdb.txt', sep="\t", index=None)
    # Fgram atts
    fn_fgram_atts = request_data(url='https://stringdb-static.org/download/protein.info.v11.0/5518.protein.info.v11.0.txt.gz', fn_out=f"{base}/string/fgram_atts.txt.gz")
    fn_fgram_atts = unzip_tidy(fn_fgram_atts, '')
    fgram_atts_df = pd.read_csv(fn_fgram_atts, sep="\t")
    fgram_atts_df['protein_external_id']  = split_df(fgram_atts_df, 'protein_external_id', '.', 1)
    fgram_atts_df.to_csv(f'{base}/string/fusarium_gram_attributes.txt', sep="\t", index=None)

    # F pgram
    fn_pseudogram_ppi_stringdb = request_data(url='https://stringdb-static.org/download/protein.links.v11.0/101028.protein.links.v11.0.txt.gz', fn_out=f"{base}/string/Fusarium_psuedo_ppi.txt.gz")
    fn_pseudogram_ppi_stringdb = unzip_tidy(fn_pseudogram_ppi_stringdb, '')
    f_psuedo_df = pd.read_csv(fn_pseudogram_ppi_stringdb, sep=" ")
    f_psuedo_df['protein1'], f_psuedo_df['protein2'] = split_df(f_psuedo_df, 'protein1', '.', 1), split_df(f_psuedo_df, 'protein2', '.', 1)
    f_psuedo_df['protein1'], f_psuedo_df['protein2'] = split_df(f_psuedo_df, 'protein1', 'P0', 0), split_df(f_psuedo_df, 'protein2', 'P0', 0)

    f_psuedo_df.to_csv(f'{base}/string/fpsuedo_ppi_stringdb.txt', sep="\t", index=None)
    # f psuedo atts
    fn_psuedo_atts = request_data(url='https://stringdb-static.org/download/protein.info.v11.0/101028.protein.info.v11.0.txt.gz', fn_out=f"{base}/string/fpsuedo_atts.txt.gz")
    fn_psuedo_atts = unzip_tidy(fn_psuedo_atts, '')
    fn_psuedo_atts_df = pd.read_csv(fn_psuedo_atts, sep="\t")

    fn_psuedo_atts_df['protein_external_id']  = split_df(fn_psuedo_atts_df, 'protein_external_id', '.', 1)
    fn_psuedo_atts_df.to_csv(f'{base}/string/fusarium_psuedo_atts.txt', sep="\t", index=None)

def egg_nog_blast_data(base):
    plant_species_list = ['nectriaceae', 'sordariomycetes', 'sordariales', 'sordariaceae',
                          'glomerellales', 'ascomycota', 'hypocreales', 'eukaryota', 'dothideomycetes',
                          'eurotiales', 'dothideomycetidae', 'clavicipitaceae', 'pleosporales',
                          'fungi', 'chaetothyriomycetidae', 'magnaporthales', 'opisthokonta', 'hypocreaceae'
                          'onygenales', 'eurotiomycetes', 'ophiostomatales', 'alphaproteobacteria', 'arthrodermataceae'
                          'agaricomycetes_incertae_sedis', 'chaetomiaceae', 'pythiales', 'thiotrichales', 'leotiomycetes'
                          'kinetoplastida', 'rhizobiaceae', 'poales', 'fusarium']

    blast_paths = [f"{base}/BLAST/species/{x}" for x in plant_species_list]

    df_blast_lists = blast_concat(paths=blast_paths, fn="f_culmorum_out.txt")
    df_eggnog_blast = pd.concat(df_blast_lists)
    df_eggnog_blast = df_eggnog_blast.drop_duplicates()
    df_eggnog_blast.to_csv(f"{base}/BLAST/f_culmorum_mapping_uniprot.txt", sep="\t", index=None, header=True)

def fusarium_gene_pro_mapping(base):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(open(f'{(base + "/eggNog")}/fculmorumUK99vs_proteins.fa'), "fasta"))

    protein_ids = list(fasta_dict.keys())
    gene_ids = [x.replace("T", "G") for x in protein_ids]

    gene_pro_mapping = {'gene ids': gene_ids,
                        'protein ids': protein_ids
                       }
    gene_pro_mapping_df = pd.DataFrame(gene_pro_mapping, columns=['gene ids', 'protein ids'])
    gene_pro_mapping_df.to_csv(f"{base}/mapping/fusarium_culmorum_gene_protein_mapping.txt", sep="\t", index=None, header=True)

def phibase_aggregate(df, agg_col):

    updated_df = df[['Gene', 'Gene ID', 'Protein ID', 'Host species', 'Pathogen species', agg_col]]
    updated_df = updated_df.replace(r'no data found', np.nan, regex=True)
    updated_df.dropna(inplace=True)
    updated_df.drop_duplicates(inplace=True)
    updated_df = updated_df.groupby(agg_col).agg({'Host species':'first',
                                                    'Protein ID': ';'.join,
                                                    'Pathogen species': ';'.join}).reset_index()

    updated_df = df_split_col_delimiter(updated_df, 'Protein ID', ';')
    updated_df = df_split_col_delimiter(updated_df, 'Pathogen species', ';')
    updated_df.drop_duplicates(inplace=True)
    updated_df.reset_index(drop=True, inplace=True)
    fusarium_names = ['venenatum', 'pseudograminearum', 'oxysporum', 'culmorum', 'langsethiae']
    print("Finished")
    updated_df = updated_df[updated_df['Pathogen species'].str.contains('|'.join(fusarium_names))]
    return updated_df[['Protein ID', 'Host species', 'Pathogen species', agg_col]]


def phibase_mapping(base):

    phi_base_blast_raw_df = pd.read_csv(f"{base}/phibase/phibase_blast_raw.out", sep="\t", header=None)
    phi_base_blast_mapping_df = pd.read_csv(f"{base}/phibase/f_culmorum_phi_mapping.txt", sep="\t")
    phi_base_blast_mapping_df.columns = ['Gene', 'Protein ID']
    phi_fn = wget.download('https://raw.githubusercontent.com/PHI-base/data/master/releases/phi-base_current.csv', out=f"{base}/phibase/")


    phi_df = pd.read_csv(phi_fn, sep=",")
    # Remove unamed columns
    phi_df.drop(phi_df.columns[phi_df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
    col_names = phi_df.columns
    col_names_list = list(col_names)
    phenotype_names = [x for x in col_names_list if 'phenotype' in x.lower()]

    phi_series = phi_df['Pathogen species'].str.upper()
    fusarium_index = phi_series[phi_series.str.contains("FUSARIUM")].index # 26 instances of Fusarium

    fusarium_phi = phi_df[phi_df.index.isin(fusarium_index)] # fusarium specific Phi results
    updated_fusarium_df = fusarium_phi[['Gene', 'Gene ID', 'Protein ID', 'Host species', 'Pathogen species', 'Disease', 'Mutant Phenotype']]
    updated_fusarium_df.reset_index(inplace=True)
    del updated_fusarium_df['index']

    disease_df = phibase_aggregate(updated_fusarium_df, 'Disease')
    disease_df.to_csv(f"{base}/phibase/fusarium-phibase-disease.txt", sep="\t", index=None)
    phenotype_df = phibase_aggregate(updated_fusarium_df, 'Mutant Phenotype')
    phenotype_df.to_csv(f"{base}/phibase/fusarium-phibase-phenotype.txt", sep="\t", index=None)

    gene_mapping_fusarium_phi_df = fusarium_phi[['Gene', 'Gene ID', 'Protein ID']]
    gene_mapping_fusarium_phi_df.columns = ['Gene name', 'Gene ID', 'Protein ID']

    phi_base_blast_mapping_df = phi_base_blast_mapping_df.drop_duplicates(subset='Gene', keep='first')


    gene_name_phibase_merged = pd.merge(phi_base_blast_mapping_df, gene_mapping_fusarium_phi_df, on="Protein ID", how='inner')
    gene_name_phibase_merged = gene_name_phibase_merged[['Gene', 'Gene name']]
    gene_name_phibase_merged = gene_name_phibase_merged.drop_duplicates(subset='Gene name', keep='first')

    updated_gene_names = []
    for name in gene_name_phibase_merged['Gene name']:
        if name.startswith("("):
            updated_gene_names.append(name.replace("(", "").replace(")", ""))
        elif "(" in name:
            updated_gene_names.append(name.split('(', 1)[0])
        else:
            updated_gene_names.append(name)

    gene_name_phibase_merged['Gene name'] = updated_gene_names
    gene_name_phibase_merged['Gene'] = gene_name_phibase_merged['Gene'].str.replace("T", "G")
    gene_name_phibase_merged.to_csv(f"{base}/phibase/fusarium-phi-gene-mapping.txt", sep="\t", index=None)

    phi_base_blast_raw_df = phi_base_blast_raw_df.drop_duplicates(subset=0, keep='first')
    phi_base_blast_raw_df.to_csv(f"{base}/phibase/phibase-blast-filtered.txt", sep="\t", index=None)


def mutant_names_fcul(base):
    fusarium_mutant_db_df = pd.read_csv(f"{base}/misc/fusarium_mutant_db.tsv", sep="\t")

    f_mutant_db_updated_Df = fusarium_mutant_db_df[['Ensembl (CS3005)', 'Gene name', 'FungiDB.1', 'Functional category of deleted gene', 'Type of mutant', 'PHI-base']]
    f_mutant_db_updated_Df['FungiDB.1'].dropna(inplace=True)
    f_mutant_db_updated_Df['Functional category of deleted gene'] = f_mutant_db_updated_Df['Functional category of deleted gene'].replace(np.nan, "Not recorded")
    f_mutant_db_updated_Df['Type of mutant'] = f_mutant_db_updated_Df['Type of mutant'].replace(np.nan, "Not recorded")

    gene_to_transcript_fg = f_mutant_db_updated_Df['FungiDB.1'].str.split("_").str[1].str.replace("G", "T")
    f_mutant_db_df_gene_list = f_mutant_db_updated_Df['FungiDB.1'].str.split("_").str[0] + "_" + gene_to_transcript_fg
    del f_mutant_db_updated_Df['FungiDB.1']
    f_mutant_db_updated_Df['FungiDB.1'] = f_mutant_db_df_gene_list


    f_mutant_db_updated_Df.to_csv(f"{base}/misc/fg_gene_names.txt", sep="\t", index=None)

def blast_2_go(base):
    blast2go_df = pd.read_csv(f"{base}/blast2go/fcul_blast2go.tsv", sep="\t", header=None)
    del blast2go_df[1]
    del blast2go_df[2]
    blast2go_df[3] = blast2go_df[3].replace("---NA---", np.nan)
    # Drop all the BLAST relationships as we already have them
    blast2go_df = blast2go_df[~blast2go_df[3].str.contains("hypothetical protein", na=False)]
    blast2go_df[3] = blast2go_df[3].replace("unamed protein product", np.nan)
    blast2go_df.dropna(inplace=True)
    blast2go_df[3] = blast2go_df[3].str.replace("Fusarium oxysporum", "F.oxy")
    blast2go_df[3] = blast2go_df[3].str.replace("fungus", "")
    blast2go_df[3] = blast2go_df[3].str.replace("  ", " ")
    #blast2go_df[3] = blast2go_df[3].str.replace("hypothetical protein ", "")
    blast2go_df.columns = ['protein', 'protein annotation']
    blast2go_df.to_csv(f"{base}/blast2go/fcul_blast2go.txt", sep="\t", index=None)
    ## updating ascomyata blast data
    asc_blast_df = pd.read_csv(f"{base}/BLAST/f_ascomycota_blast.txt", sep="\t", header=None)
    asc_blast_df = asc_blast_df.groupby(0).head(10)
    asc_blast_df.to_csv(f"{base}/f_ascomycota_blast_updated.txt", sep="\t", index=None)

def neurosporta_gene_names(base):
    # This method requires BLAST data from nuerospora linking neurospora IDs with Fusarium culmorum ID's and then getting respective gene names based on gene IDs
    # Alternatively can use OMA data - OMA link --> gene name

    # BioMart manually obtained from URL http://fungi.ensembl.org/biomart/martview/6026ccddbdebe1ecc42a394714681a77
    # neurosport_mart_fn = wget.download("http://fungi.ensembl.org/biomart/martview?VIRTUALSCHEMANAME=fungi_mart&ATTRIBUTES=ncrassa_eg_gene.default.feature_page.ensembl_gene_id|ncrassa_eg_gene.default.feature_page.ensembl_transcript_id|ncrassa_eg_gene.default.feature_page.external_gene_name&FILTERS=&VISIBLEPANEL=resultspanel", out=f"{base}/mapping/nuerospora_mart.txt")
    # Using OMA output - must be imported.
    oma_ncrassa_df = pd.read_csv(f"{base}/OMA/FculmorumvsNcrassaEnsembl-mapping.txt", sep="\t", header=None)
    # tidy up
    del oma_ncrassa_df[0]
    oma_ncrassa_df[1] = oma_ncrassa_df[1].str.split(":").str[1]
    oma_ncrassa_df[2] = oma_ncrassa_df[2].str.split(":").str[1]
    oma_ncrassa_df[2] = oma_ncrassa_df[2].str.split("_").str[0]
    oma_ncrassa_df.columns = ['Fusarium Protein ID', 'Nuerospora Protein ID']

    nuerospora_biomart_df = pd.read_csv(f"{base}/biomart/ncrassa-biomart.txt", sep="\t")
    nuerospora_biomart_df.dropna(inplace=True)
    nuerospora_biomart_df.columns = ['Gene ID', 'Nuerospora Protein ID', 'Gene name']
    del nuerospora_biomart_df['Gene ID']

    nuerospora_merged_df = pd.merge(oma_ncrassa_df,nuerospora_biomart_df, on="Nuerospora Protein ID", how='inner')
    nuerospora_merged_df.dropna(inplace=True)
    # Gene names only
    del nuerospora_merged_df['Nuerospora Protein ID']
    nuerospora_merged_df['Fusarium Protein ID'] = nuerospora_merged_df['Fusarium Protein ID'].str.replace("T", "G")
    nuerospora_merged_df = nuerospora_merged_df.drop_duplicates(subset='Gene name', keep="first")
    nuerospora_merged_df.to_csv(f"{base}/OMA/nuerospora-fcul-gene-mapping.txt", sep="\t", index=None)

def fgram_gene_names(base):

    oma_fgram_df = pd.read_csv(f"{base}/OMA/FculmorumvsFgramEnsembl-mapping.txt", sep="\t", header=None)
    # tidy up
    del oma_fgram_df[0]
    oma_fgram_df[1] = oma_fgram_df[1].str.split(":").str[1]
    oma_fgram_df[2] = oma_fgram_df[2].str.split(":").str[1]
    oma_fgram_df.columns = ['Fusarium Protein ID', 'rresv5']

    frgram_mapping_df = pd.read_csv(f"{base}/mapping/fgraminearumalias.tsv", sep="\t")

    gene_mapping_df = pd.merge(oma_fgram_df, frgram_mapping_df, on="rresv5", how="inner" )
    del gene_mapping_df['rresv5']
    gene_mapping_df.drop_duplicates(inplace=True)
    gene_mapping_df.dropna(inplace=True)
    gene_mapping_df['Fusarium Protein ID'] = gene_mapping_df['Fusarium Protein ID'].str.replace("T", "G")
    gene_mapping_df.to_csv(f"{base}/OMA/fgram-fcul-gene-name-mapping.txt", sep="\t", index=None)


#base = '/home/joseph/data'
parser = OptionParser()
parser.add_option("-b", "--bdir", type="string",
                  help="Base directory for where all files will be written to or prestored.",
                  dest="base")
parser.add_option("-e", "--ensembl", type="string",
                  help="Boolean to determine if you wish to download Ensembl data or not, t or f", dest="ensembl")
parser.add_option("-egg", "--eggnog", type="string",
                  help="Boolean to determine if you wish to download eggNog data or not, t or f", dest="egg")
parser.add_option("-m", "--mapping", type="string",
                  help="Boolean to determine if you wish to download mapping data or not, t or f", dest="map")
parser.add_option("-n", "--names", type="string",
                  help="Boolean to determine if you wish to download additional gene-name data or not, t or f", dest="name")
parser.add_option("-p", "--phi", type="string",
                  help="Boolean to determine if you wish to download phibase data or not, t or f", dest="phi")
parser.add_option("-b2g", "--blast2go", type="string",
                  help="Boolean to determine if you wish to download blast2go data or not, t or f", dest="b2g")
parser.add_option("-str", "--string", type="string",
                  help="Boolean to determine if you wish to download string data or not, t or f", dest="b2g")
options, arguments = parser.parse_args()
if options.base:
    print(f"Base directory given as {options.base}")
    base = options.base


file_dirs = ['uniprot', 'BLAST', 'cyc', 'InterPro', 'eggNog', 'mapping', 'ensembl', 'agdb', 'string', 'OMA', 'biomart']
[mkfile(f'{base}/{dir}') for dir in file_dirs] # Make the folders

if options.ensembl:
    ensembl_bool = options.ensembl
    if ensembl_bool.upper() == "TRUE" or ensembl_bool.upper() == "T":
        # Ensembl data & file names
        f_cul_df = blast_extract(path = f"{base}/BLAST/", fn = "results_f_culmorum.out", True)
        ascomycota_df = blast_extract(path = f"{base}/BLAST/", fn = "all_uniprot_f_culmorum.out", True)
        f_cul_df.to_csv(f'{base}/BLAST/f_culmorum_phi_mapping.txt', sep="\t", index=None, header=True)
        ascomycota_df.to_csv(f'{base}/BLAST/f_culmorum_ascomycota_mapping.txt', sep="\t", index=None, header=True)

if options.egg:
    egg_bool = options.egg
    if egg_bool.upper() == "TRUE" or egg_bool.upper() == "T":
        # eggNog specific data
        print("Fetching eggNog data")
        othor_taxid = egg_nog(path=(base + "/eggNog"), fn_egg="egg_nog_fusarium_filtered.tsv", fn_fasta="fculmorumUK99vs_proteins.fa")
        othor_taxid['fusarium'] = 5506
        # Grab the UniProt data
        uniprot_ftp(othor_taxid, species)
        print("Fetching data to perform BLAST with later...\n")
        egg_nog_blast_data(base)
        print("Finished fetching UniProt data for species as defined by eggNOG\n")

if options.map:
    map_bool = options.map
    if map_bool.upper() == "TRUE" or map_bool.upper() == "T":
        # Map gene ID to protein ID - for use with FASTA/GFF3 parser
        print("Fetching mapping data for Fusarium Culmorum\n")
        fusarium_gene_pro_mapping(base)

if options.name:
    name_bool = options.name
    if name_bool.upper() == "TRUE" or name_bool.upper() == "T":
        # Fetching additional names from independent dataset
        print("Fetching additional names from independent dataset and sortings them")
        mutant_names_fcul(base)
        print("Finished!\n")

if options.phi:
    phi_bool = options.phi
    if phi_bool.upper() == "TRUE" or phi_bool.upper() == "T":
        # Fetch the phibase data to add to the KG
        print("Fetching PhiBase data")
        phibase_mapping(base)
        print("Finished downloading PhiBase data")

# if options.b2g:
#     b2g_bool = options.b2g
#     if b2g_bool.upper() == "TRUE" or b2g_bool.upper() == "T":
        # BLAST2GO - Not used
        #blast_2_go(base)

if options.string:
    string_bool = options.string
    if string_bool.upper() == "TRUE" or string_bool.upper() == "T":
        string_ppi_data(base)
