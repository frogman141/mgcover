import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def proc_annotate_id(df_row):
    # this is a method to be applied to the annotated table dataframe
    _id = df_row['sequence id']

    try:
        annotate_ids = _id.split("|")[1]
        return annotate_ids
    except:
        print("skipping a row...")


def fetch_annotate_ids(annotate_table=None):
    print("fetch seq ids from annotated table...")

    annotate_table['sequence id'] = annotate_table.apply(proc_annotate_id, axis=1)

    return annotate_table['sequence id'].values


def fetch_peptide_cluster_ids(peptide_table=None, annotate_seq_ids=None):
    print ("fetch seq ids from peptide cluster data...")

    seq_id_found = []
    seq_id_not = []

    for _id in annotate_seq_ids:

        if _id in peptide_table['annotate_id'].values: # this may be extremely inefficient
            row = peptide_table.loc[peptide_table['annotate_id'] == _id]
            seq_id_found.append(_id)

            for seq_id in row['clusters_seq_ids'].values:
                seq_id_found.append(seq_id)

        else:
            seq_id_not.append(_id)
            continue

    return seq_id_found, seq_id_not


def fetch_rRNA_cluster_ids(rRNA_table=None, seq_id_not=None, seq_id_found=None):

    print("fetch seq ids from rRNA cluster data...")

    no_cluster_ids = []

    for _id in seq_id_not:

        if _id in rRNA_table['annotate_id'].values:  # this may be extremely inefficient
            row = rRNA_table.loc[rRNA_table['annotate_id'] == _id]
            seq_id_found.append(_id)

            for seq_id in row['clusters_seq_ids'].values:
                seq_id_found.append(seq_id)

        else:
            no_cluster_ids.append(_id)
            continue

    return seq_id_found, no_cluster_ids


def process_seq_ids(seq_ids=None):
    print("cleaning sequence ids:removing fragGeneScan encoding...")
    process_ids = []
    counter = 0

    for _id in seq_ids:

        try:
            fasta_formatted = _id.split("_")[0]
            process_ids.append(fasta_formatted)

        except:
            counter = counter + 1
            continue

    return process_ids


def extract_ngs_read(seq_ids=None, ngs_reads=None, taxa=None, cache_fp=None):

    print ("extracting NGS reads:post Quality Control...")

    # for each ID found in the processed ngs req data add to a fasta file
    # output of this file is the output.

    taxa_seqs = []

    for n in range(len(seq_ids)):
        seq_id = seq_ids[n]
        ngs_seq_ids = ngs_reads['seq_id'].values

        if seq_id in ngs_seq_ids:
            # creating a Sequence Record object that'll be used to create a fasta file of all of the
            # sequence that have been extracted based off of taxon annotation.

            taxa_seqs.append(
                SeqRecord(
                    ngs_reads['seq'][n],
                    id=seq_id
                    )
                )

        else:
            print ("Didn't find sequence id...")

    extract_fp = cache_fp + "extract-" + str(time.time()) + "-" + taxa + ".fa"
    SeqIO.write(taxa_seqs, extract_fp, "fasta")

    return extract_fp


def extract(annotation=None, peptide_cluster=None, rRNA_cluster=None, ngs_reads=None, taxa=None, cache=None):

    # main function of the extract module. This is used as a higher level api for extracting NGS reads from annotated tables.
    annotate_seq_ids = fetch_annotate_ids(annotate_table=annotation)
    seq_id_found, seq_id_not = fetch_peptide_cluster_ids(peptide_table=peptide_cluster,
                                                         annotate_seq_ids=annotate_seq_ids)

    # currently this function is doing nothing?
    seq_ids, no_cluster_ids = fetch_rRNA_cluster_ids(rRNA_table=rRNA_cluster, seq_id_found=seq_id_found,
                                                     seq_id_not=seq_id_not)
    seq_ids.extend(no_cluster_ids)

    process_ids = process_seq_ids(seq_ids=seq_ids)
    extract_path = extract_ngs_read(seq_ids=process_ids, ngs_reads=ngs_reads, taxa=taxa, cache_fp=cache)

    return extract_path


if __name__ == '__main__':
    from mgrast import mgrast

    print("Testing the extract module...")

    root_url = "http://api.metagenomics.anl.gov/1"
    key = "Cu3qv9ftPRaVVUYExeguMiUjR"
    id = "mgm4755753.3"
    source = "RefSeq"
    taxa = "Treponema"
    level = 'genus'
    cache = './cache/'

    annotation, peptide_cluster, rRNA_cluster, ngs_reads = mgrast(root_url=root_url, key=key, id=id,
                                                                                    source=source, taxa=taxa,
                                                                                    level=level, cache=cache)

    extract_path = extract(annotation=annotation, peptide_cluster=peptide_cluster, rRNA_cluster=rRNA_cluster, taxa=taxa, ngs_reads=ngs_reads, cache=cache)

    print (extract_path)