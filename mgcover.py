import os
import subprocess
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from urllib import urlencode
from urllib2 import urlopen, Request

"""
    MG-Cover:

    MG-Cover is a taxa specific binning to coverage calculation pipeline and is a complimentary metagenomics pipeline to 
    MG-RAST. Occasionally you'll see commits similar to this one inside the source code. These comments are meant to provide
    any individual technical explanations about the inner workings/implementation of MG-Cover.

    Thank you for being interested in my source code. I hope these commits offers some assistance.

    Sincerely,

    Alexander Baker
"""

class Mgrast(object):
    
    """
        Overview: The Mgrast class is responsible for containing all of the methods related to interacting with the MG-RAST
        API. The main method of this class is the fetch_annotation() method. And the classes outputs are the annotated table of 
        a specific Taxa and the names of the Taxa. 

        Inputs:
            - key: Web key of MG-RAST. In order for MG-Cover to function users must provided their MG-RASTs users account web key. This
            web key (also known as an API key) is what allows MG-Cover to extract annotated tables from MG-RAST.

            - id: MG-RAST run id. When you run sequences through MG-RAST an id is generated. This id is needed by MG-Cover to in order to 
            download the annotated tables of a specific MG-RAST job. 

            - taxa: The specific genus or species an individual wishes to extract

            - cache: an internal filepath to the programs cache (pending refactoring)
    """


    def __init__(self, key, id, taxa, cache):
        self.key = key
        self.id = id 
        self.taxa = taxa
        self.cache = cache
        self.root_url = "http://api.metagenomics.anl.gov/1"
    
    def error(self, msg):
        raise Exception(msg)

    def stdout_from_api(self, url=None):
        """
            Overview: stdout_from_api is responsible for sending API requests for data from the MG-RAST api. Upon recieving
            the data it's returned to calling function as a string data stream.

            Inputs:
                - url: a customily formatted URL for the MG-RAST API

            Ouputs:
                - res.read(): this function is converting the API response into a stream of text data.
        """
        header = {'Accept': 'text/plain'}
        header['Auth'] = self.key

        try:
            req = Request(url, headers=header)
            res = urlopen(req)
        except:
            raise

        if not res:
            print("ERROR: No Results where returned...")
        else:
            return res.read()

    def taxa_parser(self):
        """
            Overview: taxa_parser is responsible for parsing the inputted taxa's from the command line --taxa flag. 
            This parser utilizes : to separate different taxa's.

            Outputs:
                - a list of taxa's that can be used to request annotated tables from MG-RAST
        """

        print "Parsing Taxa Inputs..."
    
        try:
            if len(self.taxa) >= 1:
                taxa_input = " ".join(self.taxa)

                if ":" in taxa_input:
                    return taxa_input.split(":")
            
            return [taxa_input]
        
        except:
            self.error("error occurred when parsing taxa inputs...")
    
    def write_cache(self, file_name=None, file_type=None, file_data=None):
        """
            Overview: write_cache is responsible for writting caches for the MG-RAST table.

            Inputs:
                - file_name: name of the file being used for caching

                - file_type: whether it is a text stream or a normal file.

                - file_data: the actually data being written to the file
            
            Outputs:
                - there is no directly returned object. Just the filepath being written.
        """

        if file_type == 'stream':
            mode = 'wb'
        else:
            mode = 'w'

        with open(file_name, mode) as cache_file:
            cache_file.write(file_data)

    def fetch_annotation(self, source='RefSeq'):
        """
            Overview: This is the main method in the Mgrast class. It's called in the main function and is where all of the 
            the methods are assembled into a functioning block. This method utilizes all of the class inputs. This method 
            works by first parsing the taxa's input and begin looping through them. These taxa's are used to construct a MG-RAST
            API url to download a taxa specific annotated table. 

            Outputs:
                - annotated_tables: a list of file paths that lead to the annotated tables downloaded from MG-RAST

                - taxa: a list of taxa's the user inputted (need to check if I am not going to need this)

            NOTE: CURRENTLY THE PROGRAM CAN'T FUNCTION WITHOUT A TAXA INPUT MUST FIX THIS!!!!!
        """

        taxa = self.taxa_parser()
        print "Fetching Annotation Table..."
        annotated_tables = []
        url_params = [
               ('source', source),
               ('evalue', None),
               ('identity', None),
               ('length', None),
               ('type', 'organism')
            ]

        for n in taxa:
            # primitive error handling here
            try:
                url_params.append(('filter', taxa))
                url = self.root_url + "/annotation/sequence/"+ self.id + "?" + urlencode(url_params, True)
                
                if " " in n:
                    temp = n.split(" ")
                    cleaned = "_".join(temp)
                else:
                    cleaned = n

                cache_fp = self.cache + cleaned + '_annotation_table.tsv'

                annotate_stream = self.stdout_from_api(url=url)
                self.write_cache(file_name=cache_fp, file_type='stream', file_data=annotate_stream)
                annotated_tables.append(cache_fp)
            
            except Exception as e:
                raise e

        return annotated_tables, taxa

class SeqProcessing(object):

    """
        Overview: SeqProcessing class is responsible for containing all of the methods that clean and extract the sequences and sequence ids of MG-RASTs annotated table
        and their pair-end. The main method of this class is process(). It's responsible for assemblying all of the classes methods into on excuetable block of code. 
        The output of SeqProcessing are two fasta files of pair end sequences and a string to query for the reference genome.

        Inputs:
            - annotated_tables: filepath to an annotated table

            - cache: filepath the the cache location

            - pair_end: the pair end sequences. (file must be located in the current directory)

            - mgrast_id: the id of the MG-RAST run
    """

    def __init__(self, annotated_tables, cache, pair_end, mgrast_id):
        self.tables = annotated_tables
        self.cache = cache + mgrast_id
        self.mgrast_id = mgrast_id
        self.pair_end = pair_end

    def error(self, msg):
        print "Class: SeqProcessing"
        raise Exception(msg)

    def clean_ids(self, ids=None):
        """
            Overview: This method is responsible for cleaning the ids extracted froms MG-RAST annotated sequence tables. 
            This is achieved by spliting the | and _ characters used to append information onto FASTA seq ids. These cleaned
            ids are then appended to the cleaned list and returned to the process method to be used for extracting pair end sequences.

            Input:
                - ids: ids extracted from the MG-RAST annotated table
            
            Output:
                - cleaned: a python list contain cleaned FASTA seq ids from MG-RASTs annotated table.
        """
        cleaned = []
        for id in ids:
            if "|" in id:
                pre = id.split("|")[1]
                fasta_id = pre.split("_")[0]
                cleaned.append(fasta_id)
            
        return cleaned

    def extract_pair_ends(self, ids=None):
        """
            Overview: This method is reponsible for extracting the pair end FASTQ file provided via the command line.
            This extraction utilizes the clean sequence ids of the annotated tables from MG-RAST. Because pair end sequences
            use identical sequence ids we can extract sequences based off of seq ids.

            Input:
                - ids: cleaned sequence ids from MG-RASTs annotated table
            
            ouput:
                - the filepath location of the pair end FASTA file created. 
        """
        pair_end = SeqIO.index(self.pair_end, "fastq")
        pe_fileName = self.cache + self.pair_end.split(".")[0] + ".fasta"
        write_fasta = []

        for n in ids:
                write_fasta.append(pair_end[n])

        SeqIO.write(write_fasta, pe_fileName, "fasta")

        return pe_fileName

    def extract_ngs_read(self, seq_ids=None, ngs_reads=None):
        """
            Overview: This method is reponsible for creating a FASTA file for all of the sequences from the MG-RASTs annotated table being used.
            This is done by utilizing the pandas library to return a number array of the DNA sequence and sequence id from the annotated table. 
            Then this method loops through the sequence ids and creates a list of SeqRecords that are then used to create a FASTA file. These SeqRecords
            are created by BioPython library. 

            Input:
                - seq_ids: cleaned sequence ids from MG-RASTs annotated table

                - ngs_reads: pandas dataframe of the MG-RAST annotated table 
            
            ouput:
                - the filepath location of the FASTA file created. 
        """

        taxa_seqs = []
        fileName = self.cache + ".fasta"
        sequence = ngs_reads['dna sequence'].values
        ids_check = ngs_reads['sequence id'].values
        
        for indx, value in enumerate(seq_ids):

            if value in ids_check[indx]:
                taxa_seqs.append(SeqRecord(
                    Seq(sequence[indx]),
                    name=value,
                    id=value,
                    description=self.mgrast_id,
                ))

        SeqIO.write(taxa_seqs, fileName, "fasta")

        return fileName
    
    def ref_genome_search(self, annotations=None):
        """
            Overview: This method is responsible for extracting the MG-RAST assigned taxa. This assigned taxa is then utilized by
            MG-Cover to query the complete genome from NCBI. This is acheived by spliting the taxa's identified by MG-RAST. MG-RAST separates
            identified taxa's by a ';'. This method first splits the identified taxa's and then grabs the whole scientific name of the taxa.

            Input:
                - annotations: this is the list of annotations of the first row of annotated sequences by MG-RAST
            
            Output:
                - ref: a string containing the reference genome.
        """

        anno = annotations[0].split(";")
        return anno

    def processing(self):
        """
            Overview: This method is responsible for assembling all of SeqProcessing class into a functional block. This functional block first 
            processes and cleans the annotated sequence ids of MG-RAST. Then extracts the pair end sequences. Finally, the name of the reference
            genome is extracted from MG-RAST. 

            Ouput:
                - pair_end: filepath to pair end sequences of FASTQ file MG-RAST annotated

                - extract_end: filepath to the FASTA created for the sequences found in the annotated sequences of MG-RAST.

                - reference: a string containing the scientific name of the species MG-RAST identified the sequence 
        """
        
        # try:
        df = pd.read_csv(self.tables, sep='\t')
        annotated_ids = df["sequence id"].values

        ids_cleaned = self.clean_ids(ids=annotated_ids)
        pair_end = self.extract_pair_ends(ids=ids_cleaned)
        # we need to clean the pair ends with solexaQA. Currently where not trimming our adapters.

        extract_end = self.extract_ngs_read(seq_ids=ids_cleaned, ngs_reads=df)

        reference = self.ref_genome_search(annotations=df["semicolon separated list of annotations"].values)

        return pair_end, extract_end, reference




class CoverageAnalysis(object):

    """
        Overview: The CoverageAnalysis class is responsible for extracting a reference genome and conducting a 
        Coverage Anaylsis of the annotated pair end sequences extracted from MG-RAST. The main method of this class
        is the analysis(). The results of this class are the files outputed during process of calculating coverage 
        depth.

        Inputs:
            - extract: FASTA file created from the extracted sequeces of MG-RAST annotated table

            - pair: FASTA file created from the pair end sequences.

            - email: user's email to inform NCBI API of who the user is

            - taxa: user entered taxa

            - search_term: the search_term utilized the query the NCBI API Nucleotide database

            - ref_genome: the genomes MG-RAST assigned to be used to search for the reference genome

            - cache: filepath the location where all files are outputted.
    """

    def __init__(self, extract, pair, email, taxa, search_term, ref_genome, cache):
        self.r1 = extract
        self.r2 = pair
        self.sai_one = cache + 'one.sai'
        self.sai_two = cache + 'two.sai'
        self.sam = cache + 'running.sam'
        self.bam = cache + 'running_bam.bam'
        self.bam_sorted = cache + 'running.sorted.bam'

        
        self.email = email
        self.taxa = taxa
        self.search_term = search_term
        self.refs = ref_genome
        self.cache = cache
        self.refs_fp = cache
        self.file_coverage = search_term[0]
    

    def check_refs(self):
        """
            (CLASS METHOD UNDER DEVELOPMENT)

            Overivew: This method is responsible for searching the current directory for the file that contains the 
            reference genome. If this file doesn't exist or you did not enter will return boolean true to search for a
            reference genome from the NCBI Entrez Databases.abs

            Output: 
                - a boolean values that tells mg-cover whether to connect to the API or not.
        """
        if self.refs:
            print "searching file"
            return True
        else:
            return True


    def ref_genome(self):
        """
            Overview: This method is responsible for searching for reference genomes that's at least a million bp 
            in length. This is done by interacting with the NCBI Entrez API. We use the search terms extracted from
            the MG-RAST annotated table to search the NCBI Nucleotide database.

            Inputs:
                - self.search_term: this class method contains the annotated taxa's MG-RAST assigned to the annotated 
                sequences being extracted.

            Output:
                - choices: a python dictionary that contains the chosen reference genome.
        """

        Entrez.email = self.email

        for term in self.search_term:
            search_res = Entrez.esearch(db='nucleotide', term=term)
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
            
            if len(choices.keys()) < 0:
                break

        return choices

    def select_ref(self, genome_choices):

        """ 
            Overview: This method is responsible for providing the user with a list to choice from the Command Line.
            These choices represent the reference genomes that were retrieved from the NCBI Entrez API. The chosen 
            reference is then used to calculate the coverage.

            Inputs: 
                - genome_choices: a python dictionary that contains potiential reference genomes for the user to chose from

            Outputs:
                - ref_fp: filepath the reference genome that was just downloaded.
        """

        keys = []
        while True:
            print('Reference Genomes for BWA: (enter a number)\n')

            for n, key in enumerate(genome_choices.keys()):
                place = n + 1
                keys.append(key)
                print ('{}) Description: {} || Length: {}'.format(place, key, genome_choices[key]['dp_count']))

            selected = str(input('\n> '))
    
            if selected.isdigit():
                key = keys[(int(selected) - 1)]
                _key = key.split(" ")
                self.file_coverage = "_".join(_key)
                ref_fp = self.refs_fp + self.file_coverage + ".fasta"
                SeqIO.write(genome_choices[key]['seq'], ref_fp, "fasta")
                return ref_fp
            else:
                print ("Please enter a digit...\n")
    
    def check_index(self, ref_fp):
        """
            Overview: This method is responsible for indexing BWA with the reference genome.
        """
        if not os.path.isfile(ref_fp + '.sa'):
            p = subprocess.Popen(["bwa", "index", ref_fp])

    def prep_align(self, ref_fp):
        """
            Overview: This method is responsible for running the BWA 1 command.
        """
        subprocess.call(["bwa", "aln", "-f", self.sai_one, ref_fp, self.r1])
        subprocess.call(["bwa", "aln", "-f", self.sai_two, ref_fp, self.r2])

    def alignment(self, ref_fp):
        """
            Overview: This method is reponsible for running the BWA sampe command.
        """
        print "running aln algorithm"
        subprocess.call(["bwa", "sampe", '-f', self.sam, ref_fp, self.sai_one, self.sai_two, self.r1, self.r2])
        
    def convert_sort(self):
        """
            Overview: This method is responsible for converting the generated sam file to bam and then to 
            sort it.
        """
        print "converting sams"
        print self.bam
        print self.bam_sorted
        subprocess.call(['samtools', 'view', '-bS', '-o', self.bam, self.sam]) 
        subprocess.call(['samtools', 'sort', '-o', self.bam_sorted, self.bam])

    def coverage_calc(self):
        """
            Overview: This method is responsible for calculating the depth of coverage via the samtools depth 
            command.
        """
        print "calculating coverage"
        coverage_file = self.cache + self.file_coverage + '.coverage'
        with open(coverage_file, 'wb') as cov:
            cov.write(subprocess.check_output(['samtools', 'depth', self.bam_sorted]))
        
        df = pd.read_csv(coverage_file, names=['genomes', 'loci', 'depth'], sep='\t')
        df.to_csv(coverage_file, index=False)

    def analysis(self):
        """
            Overview: This is the main method of the CoverageAnalysis class. Responsible for assembling all of the methods 
            from class into an excuetable block of code. First it checks for the reference genome locally then interacts with
            NCBIs API. From that it then runs the command line tools BWA and Samtools to calculate the genome coverage of the pair end
            annotated sequences.

            Output:
                - all of the files from BWA and Samtools and a coverage file that's named after the reference genome
        """
        ref_fp = self.check_refs()

        if not type(ref_fp) is str:
            choice = self.ref_genome()
            ref_fp = self.select_ref(choice)

        self.check_index(ref_fp)
        self.prep_align(ref_fp)
        self.alignment(ref_fp)
        self.convert_sort()
        self.coverage_calc()

        



def mgcover(key=None, taxa=None, id=None, pair_end=None, email=None, ref_genome=None):
    """

        Overview: This is the main function of mgcover it's responsible for creating an excutable block.

        Inputs:
            - key: users MG-RAST web key

            - taxa: names of the taxa(s) to extract. For multiple taxas separate with ":".'

            - id: MG-RAST job id for annotated sequences

            - pair_end: Email account we report to NCBIs API.

            - email: The file containing pair end sequences of the MG-RASTs annotated sequences

            - ref_genome: Fasta file containing reference genome

        Outputs:
            - a coverage file for the taxa that was specified.
    """
    cache = os.getcwd() + os.sep
    ref_genomes_dir = cache

    if id is None or pair_end is None or key is None:
        print "mgcover requires a MG-RAST job id | AND | the pair end's fastq filepath | AND | the users MG-RAST webkey..."
    
    else:
        # fetch annotated table by Taxa
        mgrast = Mgrast(key=key, id=id, taxa=taxa, cache=cache)
        annotated_tables, taxa_names = mgrast.fetch_annotation()

        # process the sequences
        for n, table in enumerate(annotated_tables):
            seqproc = SeqProcessing(annotated_tables=table, cache=cache, pair_end=pair_end, mgrast_id=id)
            pair, extract, ref = seqproc.processing()

            # Coverage Anaylsis
            coverage = CoverageAnalysis(extract=extract, pair=pair, email=email, taxa=taxa_names[n], search_term=ref, ref_genome=ref_genome, cache=cache)
            coverage.analysis()

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='mg-cover command line interface.')

    parser.add_argument('--key', type=str, help='users MG-RAST web key', default="Cu3qv9ftPRaVVUYExeguMiUjR")
    parser.add_argument('--taxa', type=str, nargs='*', help='names of the taxa(s) to extract. For multiple taxas separate with ":".', default=['Treponema'])
    parser.add_argument('--id', type=str, help='MG-RAST job id for annotated sequences')
    parser.add_argument('--email', type=str, help='Email account we report to NCBIs API.', default='bakeralex664@gmail.com')
    parser.add_argument('--pair-end', help='The file containing pair end sequences of the MG-RASTs annotated sequences.') # this is a required command currently
    parser.add_argument('--ref-genome', type=str, help='Fasta file containing reference genome')

    args = vars(parser.parse_args())
    mgcover(**args)
