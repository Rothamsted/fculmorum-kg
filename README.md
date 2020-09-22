
## General
This repository contains the ETL/scripts used to obtain data for the _F. culmorum_ KnetMiner knowledge graph. 

**N.b.** Some data is created or obtained externally, as mentioned in the manuscript; this includes eggNog data, BLAST, and OMA data. These datasets are also processed in the ```f_culmorum_etl.py``` script present in this repository.

# Dependencies

[Python 3.7+](https://www.python.org/downloads/) is a requirement for both scripts, as is [R](https://www.r-project.org/) for rpy2.

Additional Python dependencies include the following:

1. [Numpy](https://pypi.org/project/numpy/)
2. [requests](https://pypi.org/project/requests/)
3. [wget](https://pypi.org/project/wget/)
4. [gzip](https://pypi.org/project/gzip/)
5. [Pandas](https://pypi.org/project/pandas/)
6. [BioPython](https://pypi.org/project/biopython/)
7. [MySQL-python](https://pypi.org/project/mysql/)
8. [urllib3](https://pypi.org/project/urllib3/)
9. [rpy2](https://pypi.org/project/rpy2/)

# How to use the scripts

For the [biomart script](https://github.com/Rothamsted/fculmorum-kg/blob/master/fusairum_biomart.py), simply provide the directory you wish to have your data exported to. Perform the following command (append to the directory arg as you wish):

``` python fusarium_biomart.py /home/ ```

For the ```f_culmorum_etl.py``` [script](https://github.com/Rothamsted/fculmorum-kg/blob/master/f_culmorum_etl.py), which is the main ETL script, a few booleans are used to have a choice between what data you obtain. 

Please change the input file names accordingly to any updated data, directly in the script, or write a config file for it. Some of the boolean flags are grouped due to data dependencies for certain data types to be obtained/produced, i.e. UniProt data required for mapping to eggNog data. 

1. You must add a base directory with the ```-b``` or ```--bdir``` flag, which is where all your files will be written to, or respectively stored for BLAST/eggNog/OMA outputs derived externally. ^^ 
2. If you wish to download Ensembl specific data, use the ```-e``` or ```---ensembl``` flag. **NOTE:** this requires you to have your BLAST data present in the BLAST folder within your base directory, as this data is further mapped. You will need to name your data accordingly, to what's present in the script - for the _F. culmorum_ blast output, it should be named ```results_f_culmorum.out```, the _Ascomycota_ BLAST output should be named ```all_uniprot_f_culmorum.out```. The resultant outputs given are the ```f_culmorum_phi_mapping.txt``` and ```f_fulculmorum_ascomycota_mapping.txt``` outputs, containing mappings between the BLAST data and Ensembl data.
3. If you wish to specifically transformed any eggNog data (as described in the manuscript), use the ```-egg``` or ```--eggnog``` flag. You will need the eggnog output, named as ```egg_nog_fusarium_filtered.tsv```, which should be present in the eggNog folder. You will also need the base fasta file, which in this case is ```fculmorumUK99vs_proteins.fa```. You will retrieve data back in the BLAST folder, containg mappings between eggNog & BLAST proteins with gene names. Additionally, the uniprot data used will be downloaded into the uniprot directory.
4. For mapping only data use to ```-m``` or ```--mapping``` flag. This will obtain data which maps identifiers, be it gene or protein identifiers, to external identifiers. You will need to have performed the eggNog data transformation, first. 
5. For additional gene name data, use the ```-n``` or ```--names``` flag. This does require the fusarium_mutant_db.tsv curated file from RRes, which must be present in the misc folder within the base directory. The output will be called fg_gene_names.txt in the misc folder.
6. For mapping BLAST data (_F.culmorum_ to [PhiBase fasta](http://www.phi-base.org/)) to corresponding proteins in the core [PhiBase database ](https://raw.githubusercontent.com/PHI-base/data/master/releases/phi-base_current.csv), so that phenotypes and diseases can be identified, use the ```-p``` or ```--phi``` flag. Note that the input files must be present in the phibase folder within the base directory, with the blast data being named as ```phibase_blast_raw.out``` and the mapping between phibase & _F.culmorum_ named as ```f_culmorum_phi_mapping.txt```. The resultant files will be in the phibase folder,  which includes ```fusarium-phi-gene-mapping.txt``` & ```phibase-blast-filtered.txt```

^^ Note that the folder structure is indeed created by the ETL script, but you may wish to create it prior to use so dependent files can be placed in their respective folders (as outlined in the above options avaiable). 

This includes:

 'uniprot', 'BLAST', 'cyc', 'InterPro', 'eggNog', 'mapping', 'ensembl', 'agdb', 'string', 'OMA', 'biomart'
