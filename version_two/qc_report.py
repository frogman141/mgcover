from Bio import SeqIO

def parse_fasta(fasta_file=None):

    # This file is responsible for parsing fasta files

    prep_seq = {"seq_id": [], "seq": []}

    for seq in SeqIO.parse(fasta_file, "fasta"):
        prep_seq["seq_id"].append(seq.id)
        prep_seq["seq"].append(seq.seq)

    return prep_seq

if __name__ == '__main__':

    original = './job_mgm4755754.3/deruplication_one.fasta'
    solexaQA = './job_mgm4755754.3/preprocessing_two.fasta'
    drisee = './job_mgm4755754.3/deruplication_two.fasta'
    bowtie = './job_mgm4755754.3/screening_one.fasta'

    files = [solexaQA, drisee, bowtie, original]

    for n in files:

        data = parse_fasta(n)

        print ('%d: %s' % (len(data['seq']), n))