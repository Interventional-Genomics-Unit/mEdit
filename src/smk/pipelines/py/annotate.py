# Native Modules
from subprocess import Popen, PIPE
# Installed Modules
from Bio.Seq import Seq

########
###      Transcript Class is used to compute/retrieve transcript information based on a list of coordinates
###      1: Load list of genomic coordinates and ncbiRefSeq.bed file path
###      2: Call upon any of these coordinates using the NCBI transcript ID or coordinates
###         to retreive  gene names, features, reading frame, feature positions and transcript sequence
###         Created by T.Hudson
########

class Transcript:
    tx_lib = {}
    coord2tid = {}
    labels = ['chrom', 'txStart', 'txEnd', 'strand', 'tid', 'eid', 'name',
              'cdsStart', 'cdsEnd', 'exonStarts', 'exonEnds', 'exonFrames']

    def __init__(self,entry):
        '''
        :param entry: NCBI Transcript Entry
        '''

        self.entry = entry
        for k in ['txStart', 'txEnd', 'cdsStart', 'cdsEnd']:
            if k in self.entry.keys():
                self.entry[k] = int(self.entry[k])

        self.overlapping_transcripts = []

        ##output
        # mapping transcript features (relative to transcript start)
        if 'exonStarts' in self.entry.keys():
            self.exons = self.get_exons()
            self.utrs = self.get_utrs()
        self.tx_len = self.entry['txEnd'] - self.entry['txStart']
        self.flanking = [(-50, 0), (self.tx_len, self.tx_len + 50)]

        self.feature = None
        self.rf = None
        self.cdsseq = None
        self.txseq = None

    @classmethod
    def load_transcripts(cls, annote_path, snvcoords):
        '''
        Must be initiated before annotating!
        Loads coordinates into Transcript Class
        '''
        bed_data = ""
        bedtools_out =""

        # reset dict
        cls.coord2tid = {}
        cls.tx_lib = {}

        # with open(temp_bedfname, 'w') as tempbed:
        for coord in snvcoords:
            coord_field = coord.split(':')
            chrom = coord_field[0]
            start = coord_field[1].split('-')[0]
            end = coord_field[1].split('-')[1]
            line = "\t".join([chrom, start, end])
            bed_data += line + "\n"


        # Sort bedfile
        sort_cmd = ['bedtools', 'sort', '-i', "stdin"]

        sorted_bed = Popen(
            sort_cmd,
            stdin=PIPE,
            stdout=PIPE,
            text=True
        )
        sorted_output, sort_error = sorted_bed.communicate(input=bed_data)
        sorted_bed.wait()

        #find closest ORFS
        # Check for errors in the sorting step
        if sort_error:
            print("Error during sorting:", sort_error)
        else:
            # Pipe the sorted BED data directly to `bedtools window`
            window_cmd = ["bedtools", "window", "-w","50", "-a",  "stdin", "-b", annote_path]
            window_output = Popen(
                window_cmd,
                stdin=PIPE,
                stdout=PIPE,
                stderr=PIPE,
                text=True
            )

            # Run the command and capture the output
            bedtools_out, stderr = window_output.communicate(input=sorted_output)
            window_output.wait()

            if window_output.returncode == 0:
                pass
            else:
                print("Error:", stderr)

        if bedtools_out != '':

            bed_entries = bedtools_out.split('\n')[:-1]

            for bed_entry in bed_entries:
                tokens = bed_entry.split('\t')
                tokens = tokens[:-1] + tokens[-1].split('|')
                coord = f'{tokens[0]}:{tokens[1]}-{tokens[2]}'
                entry = dict(zip(cls.labels, tokens[-12:]))

                if entry['chrom'] != '.':
                    if tokens[0] == entry['chrom']:
                        if coord not in cls.coord2tid.keys():
                            cls.coord2tid[coord] = entry['tid']
                            cls.tx_lib[entry['tid']] = cls(entry)
                        # adjusting for overlaps to retrieve the most up to date transcript
                        else:
                            new_tid = entry['tid']
                            oldtid = cls.coord2tid[coord]

                            if new_tid.startswith('NR') and oldtid.startswith('NM'):
                                pass



                            elif new_tid[0:2] == oldtid[0:2]:
                                if entry['eid'] != '-' or (float(oldtid.split('_')[-1]) < float(new_tid.split('_')[-1])):
                                    cls.coord2tid[coord] = new_tid
                                    obj = cls(entry)
                                    obj.overlapping_transcripts.append(oldtid)
                                    cls.tx_lib[entry['tid']] = obj
                            else:
                                cls.coord2tid[coord] = new_tid
                                obj = cls(entry)
                                obj.overlapping_transcripts.append(oldtid)
                                cls.tx_lib[entry['tid']] = obj



    @classmethod
    def transcript(cls, coord):
        # Call upon specific coords that were loaded with bedtools
        if coord in cls.coord2tid.keys():
            tid = cls.coord2tid[coord]
            start, end = coord.split(':')[1].split('-')
            pos = int((int(start) + int(end)) /2)
            obj = cls.tx_lib[tid]
            obj.find_feature(pos)
            return obj

        else: # bedtools never annotated
        #    obj = cls('intergenic')
            return 'intergenic'

    def __dict__(self):
        return self.entry

    def get_utrs(self):
        exon_starts = self.entry['exonStarts'][:-1].split(',')
        exon_ends = self.entry['exonEnds'][:-1].split(',')
        exon_frames = self.entry['exonFrames'].replace("\n", "")[:-1].split(',')
        utrs = None

        if self.exons:
            ogexons = [(int(exon_starts[i]) - self.entry['txStart'], int(exon_ends[i]) - self.entry['txStart']) for i in range(len(exon_ends))]

            utrs = [(ogexons[0][0], self.exons[0][0] - 1),(self.exons[-1][1], ogexons[-1][1] - 1)]

            if exon_frames[-1] == '-1':  # -1 means entire exon is UTR
                utrs[1] = ogexons[-1]
            if exon_frames[0] == '-1':
                utrs[0] = ogexons[0]

        return utrs

    def get_exons(self):
        '''
        uses entry info to find relative exons positions
        '''
        exon_starts = self.entry['exonStarts'][:-1].split(',')
        exon_ends = self.entry['exonEnds'][:-1].split(',')
        exon_frames = self.entry['exonFrames'].replace("\n", "")[:-1].split(',')
        tx_start = self.entry['txStart']

        exons = [(int(exon_starts[i]) - tx_start, int(exon_ends[i]) - tx_start) for i in range(len(exon_ends))]
        for i in range(len(exon_frames)):
            if exon_frames[i] == '-1':  # -1 means entire exon is UTR
                exons = exons[1:]
                exon_starts = exon_starts[1:]
            else:
                break
        for i in range(1, len(exon_frames)):
            if exon_frames[-i] == '-1':
                exons = exons[0:len(exons) - 1]
                exon_ends = exon_ends[0:-1]
            else:
                break

        # Determine the stop and start of UTR
        if len(exons) > 0:
            exons[0] = (int(self.entry['cdsStart']) - int(exon_starts[0]) + exons[0][0], exons[0][1])
            exons[-1] = (exons[-1][0], exons[-1][1] - (int(exon_ends[-1]) - int(self.entry['cdsEnd'])))
        else:
            exons = None
        return exons


    def get_tx_seq(self,fasta_seq):
        '''
        Using a Refseq Transcript_ID, Ensembl Transcript_ID or coordinates find transcript annotations and transcript sequence
        from either a genome fasta path or given genome sequence
        '''
        self.tx_seq = fasta_seq.seq[self.entry['txStart']:self.entry['txEnd']]
        return self.tx_seq

    def get_cdsseq(self):
        '''
        uses entry info to adjust exons positions relative to transcription start
        removes utrs from exons
        translates into cds
        '''
        # Determine the stop and start of UTR
        if self.tx_seq:
            cds = Seq(''.join([str(self.tx_seq)[a:b] for a, b in self.exons]))
            if self.entry['strand'] == '-':
                cds = cds.reverse_complement()
        else:
            cds = None
        return cds


    def tx_info(self):
        return self.entry['eid'], self.entry['tid'], self.entry['name'], self.entry['strand'], self.entry['txStart']


    def find_reading_frame(self, dist_from_cds_start):
        '''
        Finds reading frame of SNV in extracted sequence
        '''
        rf = 1 if dist_from_cds_start % 3 == 2 else 2 if dist_from_cds_start % 3 == 0 else 0
        return rf

    def find_feature(self, pos):
        feature, rf = None, None
        t_snvpos = int(pos) - self.entry['txStart']

        if self.entry['tid'].startswith('NR'):
            feature = 'non-coding RNA'

        elif 'cdsStart' not in self.entry.keys():
            feature = 'Not Found'

        elif self.exons == None or t_snvpos < -50 or t_snvpos > self.tx_len + 50:
            # not in transcript - shouldn't happen or else no entry would be found
            feature = 'non-coding'

        elif t_snvpos in range(self.flanking[0][0],self.flanking[0][1]+1) or t_snvpos in range(self.flanking[1][0],self.flanking[1][1]+1):

            if t_snvpos in range(self.flanking[0][0],self.flanking[0][1]+1):
                feature = 'flanking-upstream' if self.entry['strand'] == '+' else 'flanking-downstream'
            else:
                feature = 'flanking-downstream' if self.entry['strand'] == '+' else 'flanking-downstream'


        elif t_snvpos in range(self.utrs[0][0],self.utrs[0][1]+1) or t_snvpos in range(self.utrs[1][0],self.utrs[1][1]+1):

            if t_snvpos in range(self.utrs[0][0],self.utrs[0][1]+1):
                feature = '5utr' if self.entry['strand'] == '+' else '3utr'
            else:
                feature = '3utr' if self.entry['strand'] == '+' else '5utr'


        else:# find if exon or intron
            feature = 'intron'
            exon_n = 0
            for x in self.exons:
                # if in exon find reading frame
                if t_snvpos in range(x[0], x[1] + 1):
                    exon_trans_num = exon_n + 1 if self.entry['strand'] == '+' else (len(self.exons) - exon_n)
                    feature = f'exon {str(exon_trans_num)}'
                    dist = sum([e[1] - e[0] for e in self.exons[0:exon_n]])
                    dist_from_cds_start = dist + (t_snvpos - x[0])
                    len_cds = sum([e[1] - e[0] for e in self.exons])

                    if self.entry['strand'] == '-':
                        dist_from_cds_start = (len_cds - dist_from_cds_start) + 1

                    if dist_from_cds_start < 3:
                        feature = 'start_codon'

                    if dist_from_cds_start > len_cds - 3:
                        feature = 'stop_codon'

                    rf = self.find_reading_frame(dist_from_cds_start)

                    break
                exon_n += 1

        self.feature, self.rf = feature, rf


