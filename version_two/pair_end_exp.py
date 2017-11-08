import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from urllib import urlencode
from urllib2 import urlopen, Request

def stdout_from_api(url=None, auth=None):
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
        req = Request(url, headers=header)
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

def cache_rec(cache=None):
    return pd.read_csv(cache + 'annotation_table.tsv', sep='\t')

def extract_pair_ends(fp='testing_R2.fastq', ids=None):
    pair_end = SeqIO.index(fp, "fastq")
    write_fasta = []

    for n in ids:
            write_fasta.append(pair_end[n])

    SeqIO.write(write_fasta, "taxa_extract_R2.fa", "fasta")

def clean_ids(ids=None):
    cleaned = []

    for id in ids:
        if "|" in id:
            pre = id.split("|")[1]
            fasta_id = pre.split("_")[0]
            cleaned.append(fasta_id)
        
    return cleaned


def extract_ngs_read(seq_ids=None, ngs_reads=None, taxa=None, cache_fp=None):
    taxa_seqs = []

    sequence = ngs_reads['dna sequence'].values
    ids_check = ngs_reads['sequence id'].values
    
    for indx, value in enumerate(seq_ids):

        if value in ids_check[indx]:
            taxa_seqs.append(SeqRecord(
                Seq(sequence[indx]),
                name=value,
                id=value,
                description=taxa,
            ))

    SeqIO.write(taxa_seqs, 'taxa_extract_R1.fa', "fasta")


if __name__ == '__main__':
    print("Testing the mgrast module...")
    root_url = "http://api.metagenomics.anl.gov/1"
    key = "Cu3qv9ftPRaVVUYExeguMiUjR"
    id = "mgm4763129.3"
    source = "RefSeq"
    taxa = "Treponema"
    level = 'genus'
    cache = './cache/pair_end_exp/'

    annotated_df = fetch_annotation(root_url=root_url, key=key, id=id, source=source, taxa=taxa, level=level, identity=None, evalue=None, length=None, cache=cache)
    # annotated_df = cache_rec(cache=cache)
    seqIds = clean_ids(ids=annotated_df["sequence id"].values)
    extract_pair_ends(ids=seqIds)
    extract_ngs_read(seq_ids=seqIds, ngs_reads=annotated_df, taxa=taxa, cache_fp='./')
    print ("Testing finished...")