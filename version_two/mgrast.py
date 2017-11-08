import pandas as pd
from Bio import SeqIO
from urllib.request import urlopen, Request
from urllib.parse import urlencode


def stdout_from_api(url=None, data=None, auth=None):
    """
        Overview: stdout_from_api is responsible for sending API requests for data from the MG-RAST api. Upon recieving
        the data it's returned the calling function.

        Params:
            - url: the MG-RAST REST API url that data is being requested from.

            - data: any potiential data needed for the request.

            - auth: the api_key or web_key for the MG-RAST api.

        Output:
            - file: a pandas dataframe representation of the file being returned by the api.
    """
    header = {'Accept': 'text/plain'}

    if auth:
        header['Auth'] = auth

    try:
        req = Request(url, data, headers=header)
        res = urlopen(req)
    except:
        raise

    if not res:
        print("ERROR: No Results where returned...")
    else:
        return res.read()


def write_cache(file_name=None, file_type=None, file_data=None):

    if file_type == 'stream':
        mode = 'wb'
    else:
        mode = 'w'

    with open(file_name, mode) as cache_file:
        cache_file.write(file_data)

def fetch_annotation(root_url=None, key=None, id=None, source=None, taxa=None, level=None, identity=None, evalue=None, length=None, cache=None):
    """
        Overview: This function is responsible for fetching the tabular annotation file that is the result of BLAT
        search. This search could of used RefSeq or another source of data. This is achieved through interacting with
        MG-RASTs API. This API is used to either fetch the entire annotation file or to extract specific taxons.

        Params:
            - root_url: the base url for the MG-RAST API (i.e. http://api.metagenomics.anl.gov/1)

            - key: api key (also known as web key) that mg-rast provides

            - id: MG-RAST job id. This is what allows us to extract files from a particular MG-RAST run

            - source: The data source used to annotate and extract taxa data. An example of a data source is Refseq

            - taxa: The specific taxa the user desires to be extracted. All can be answer.

            - level: The specific level in the taxa where extract sequences on (i.e. Family, Genus, etc.)

            - identity: The % of sequence identity.

            - evalue: the threshold a sequences evalue must be equal to or below.

            - cache: the file path to the cache directory.

        Output:
            - annotate_df: A pandas dataframe representation of the annotation table retrieved from MG-RAST.
    """

    print ("Fetching Annotation Table...")
    url_params = [
               ('source', source),
               ('evalue', evalue),
               ('identity', identity),
               ('length', length),
               ('type', 'organism')
            ]

    if taxa:
        url_params.append(('filter', taxa))

    url = root_url + "/annotation/sequence/"+ id + "?" + urlencode(url_params, True)
    cache_fp = cache + 'annotation_table.tsv'

    annotate_stream = stdout_from_api(url=url, auth=key)
    write_cache(file_name=cache_fp, file_type='stream', file_data=annotate_stream)

    return pd.read_csv(cache_fp, sep='\t')

def fetch_files(root_url=None, key=None, id=None, cache=None, files=None):

    # this function is responsible for fetching all of the files needed from MG-RAST
    print("Fetching all data files required...")

    dfs = {}

    for file in files:
        file_name = file[0]
        url = root_url + "/download/" + id + "?" + "file=" + file_name
        cache_fp = cache + file_name.split(".")[0] + '_table.tsv'

        file_stream = stdout_from_api(url=url, auth=key)
        write_cache(file_name=cache_fp, file_type='stream', file_data=file_stream)

        if file[1] == 'fasta':
            fasta_dict = parse_fasta(fasta_file=cache_fp)
            dfs[file_name] = pd.DataFrame.from_dict(fasta_dict)

        else:
            dfs[file_name] = pd.read_csv(cache_fp, sep='\t', names=['annotate_id', 'clusters_seq_ids', 'match_percent'])

    return dfs


def parse_fasta(fasta_file=None):

    # This file is responsible for parsing fasta files

    prep_seq = {"seq_id": [], "seq": []}

    for seq in SeqIO.parse(fasta_file, "fasta"):
        prep_seq["seq_id"].append(seq.id)
        prep_seq["seq"].append(seq.seq)

    return prep_seq


def mgrast(root_url=None, key=None, id=None, source=None, taxa=None, level=None, cache=None):

    files_needed = [['550.1', 'cluster'], ['440.1', 'cluster'], ['299.1', 'fasta']]

    annotation = fetch_annotation(root_url=root_url, key=key, id=id, source=source, taxa=taxa, level=level, cache=cache)
    dfs_of_files = fetch_files(root_url=root_url, key=key, id=id, cache=cache, files=files_needed)

    return annotation, dfs_of_files['550.1'], dfs_of_files['440.1'], dfs_of_files['299.1']

if __name__ == '__main__':
    print("Testing the mgrast module...")
    root_url = "http://api.metagenomics.anl.gov/1"
    key = "Cu3qv9ftPRaVVUYExeguMiUjR"
    id = "mgm4755753.3"
    source = "RefSeq"
    taxa = "Treponema"
    level = 'genus'
    cache = './cache/'

    annotation, peptide_cluster, rRNA_cluster, ngs_reads = mgrast(root_url=root_url, key=key, id=id, source=source, taxa=taxa, level=level, cache=cache)

    print (ngs_reads.head())
    print("mgrast module test ended...")