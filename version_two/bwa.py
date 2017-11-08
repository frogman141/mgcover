from Bio import Entrez, SeqIO
import subprocess

def ref_genome(search, email='bakeralex664@gmail.com'):
    """This is a horrible draft of this function.

    ref_genome is responsible for fetching the bwa reference genome if it's not found in our reference directory.
    """


    # setting up configs for NCBI api
    Entrez.email = email
    search_term = '"{}"[organism] AND "complete genome"[status]'.format(search)
    search_res = Entrez.esearch(db='nucleotide', term=search_term)
    res_parsed = Entrez.read(search_res)
    genome_ids = res_parsed['IdList'] # got the id's
    choices = {}
    # looping through the id's

    for genome_id in genome_ids:
        handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype='fasta')
        record = SeqIO.read(handle, "fasta")
        bp_count = len(record.seq)

        if bp_count < 1000000:
            continue

        if len(choices.keys()) - 1 < 4:
            choices[record.description] = {'seq': record , 'dp_count': bp_count, 'descriptions': record.description}
        else:
            break

    return choices

def select_ref(genome_choices):

    """This function needs to be rethought."""

    keys = []
    while True:

        print('Reference Genomes for BWA: (enter a number)\n')

        for n, key in enumerate(genome_choices.keys()):
            place = n + 1
            keys.append(key)
            print ('{}) Description: {} || Length: {}'.format(place, key, genome_choices[key]['dp_count']))

        selected = input('\n> ')
  
        if selected.isdigit():
            key = keys[(int(selected) - 1)]
            return genome_choices[key]
        else:
            print ("Please enter a digit...\n")


def run_bwa():
    print ("running bwa")


if __name__ == '__main__':
    search = "treponema pallidum"
    choice = ref_genome(search="treponema pallidum")
    reference = select_ref(choice)

    SeqIO.write(reference['seq'], './bwa_reference/%s.fa' % "treponema pallidum", "fasta")