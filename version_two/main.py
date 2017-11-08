# import pandas as pd
# from mgrast import mgrast
# from extract import extract

def result():
    print ("returning the result of mg-cover")


def checking_index(taxa=None):
    print ("checking BWA indexing records...")

    index_list = pd.read_csv("./cache/bwa_index.txt")

    # for genenome in taxa:
    #     indexed = index_list['indexed_genome'].values
    #
    #     if indexed[genenome]:


def taxa_parser(taxa):

    taxa_input = " ".join(taxa)

    if ":" in taxa_input:
        return taxa_input.split(":")
    else:
        return [taxa_input]


def main(key=None, id=None, source=None, taxa=None, level=None):

    # main function of the program entirely
    root_url = "http://api.metagenomics.anl.gov/1"
    cache = './cache/'
    extract_paths = []
    print taxa
    species_name = taxa_parser(taxa)
    print species_name
    # for species in species_name:

    #     # loop through all of the taxonomy names provided.

    #     annotation, peptide_cluster, rRNA_cluster, ngs_reads = mgrast(root_url=root_url, key=key, id=id,
    #                                                                   source=source, taxa=species,
    #                                                                   level=level, cache=cache)

    #     extract_paths.append(extract(annotation=annotation, peptide_cluster=peptide_cluster, rRNA_cluster=rRNA_cluster, ngs_reads=ngs_reads, taxa=species, cache=cache))



    # print ("Finished running mg-cover")

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='mg-cover command line interface.')

    parser.add_argument('--key', type=str, help='users MG-RAST web key', default='Cu3qv9ftPRaVVUYExeguMiUjR')
    parser.add_argument('--id', help='MG-RAST job id')
    parser.add_argument('--source', type=str, help='which database does MG-RAST utilize for annotation', default='RefSeq')
    parser.add_argument('--taxa', type=str, nargs='*', help='names of the taxa(s) to extract. For multiple taxas separate with ":".', default='Treponema')

    args = vars(parser.parse_args())
    main(**args)