import os.path

COLUMNS = ["geneName", "name", "chrom", "strand", "txStart", "txEnd",
           "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds"]

NO_EXON_COLUMNS = ["geneName", "name", "chrom", "strand", "txStart", "txEnd",
           "cdsStart", "cdsEnd", "exonCount"]

NUMERIC_COLUMNS = ["txStart", "txEnd", "cdsStart", "cdsEnd"]

NUMERIC_LIST_COLUMNS = ["exonStarts", "exonEnds"]

STRING_COLUMNS = set(COLUMNS) - set(NUMERIC_COLUMNS) - set(NUMERIC_LIST_COLUMNS)
#DELETE - ends generics
class Exon(object):
    """
    This class defines an exon inside a record
    """

    __slots__ = ["_gene", "_transcript", "_chr", "_start", "_end", "_number"]

    def __init__(self, gene, transcript, chr, start, stop, n):
        self._gene = gene
        self._transcript = transcript
        self._chr = chr
        self._start = start
        self._end = stop
        self._number = n

    @property
    def gene(self):
        return self._gene

    @property
    def transcript(self):
        return self._transcript

    @property
    def chr(self):
        return self._chr

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._end

    @property
    def number(self):
        return self._number

    @classmethod
    def fromrecord(cls, record):
        exons = []
        assert(len(record.exonStarts) == len(record.exonEnds))
        if record.strand == "+":
            for i, (s, e) in enumerate(zip(record.exonStarts,
                                           record.exonEnds)):
                exons.append(Exon(record.gene, record.transcript,
                                  record.chromosome, s, e, i+1))

        # for negative strand transcripts
        else:
            i = len(record.exonStarts)
            for (s, e) in zip(record.exonStarts, record.exonEnds):
                exons.append(Exon(record.gene, record.transcript,
                                  record.chromosome, s, e, i))
                i -= 1

        return exons


# class Transcript(object):
#     __slots__ = ["name", "gene", "chr", "start", "end", "cds_start", "cds_end", "exons", "strand"]

#     def __init__(self, name, chr, start, end, cds_start, cds_end, exons=None, gene=None, strand="+"):
#         self.name = name
#         self.gene = gene
#         self.chr = chr
#         self.start = start
#         self.end = end
#         self.cds_start = cds_start
#         self.cds_end = cds_end
#         self.exons = exons
#         self.strand = strand

#     def update_exons(self, exon):
#         if exon.start < self.start:
#             raise ValueError("Start of exon cannot be in front of start of transcript")
#         if exon.stop > self.end:
#             raise ValueError("End of exon cannot be behind end of transcript")

#         if self.exons:
#             self.exons.append(exon)
#         else:
#             self.exons = [exon]

#     @property
#     def cds_exons(self):
#         """
#         Return those exons which lie within the cds
#         Also returns those partially inside the cds
#         :return:
#         """
#         return [x for x in self.exons if x.stop >= self.cds_start and
#                 x.start <= self.cds_end]

#     @property
#     def line(self):
#         line = []
#         d = self.to_dict()
#         for nlc in NUMERIC_LIST_COLUMNS:
#             d[nlc] = ",".join(map(str, d[nlc])) + ","
#         for col in COLUMNS:
#             line += [d[col]]
#         return "\t".join(map(str, line))

#     def to_dict(self):
#         d = {}
#         d["geneName"] = self.gene.name
#         d["name"] = self.name
#         d["chrom"] = self.chr
#         d["strand"] = self.strand
#         d["txStart"] = self.start
#         d["txEnd"] = self.end
#         d["cdsStart"] = self.cds_start
#         d["cdsEnd"] = self.cds_end
#         d["exonStarts"] = [int(x.start) for x in self.exons]
#         d["exonEnds"] = [int(x.stop) for x in self.exons]
#         d["exonCount"] = len(self.exons)

#         return d


# class Gene(object):
#     __slots__ = ["name", "min_coord", "max_coord", "transcripts", "chr"]

#     def __init__(self, name, chr=None, min_coord=None, max_coord=None, transcripts=None):
#         self.name = name
#         self.min_coord = min_coord
#         self.max_coord = max_coord
#         self.transcripts = transcripts
#         self.chr = chr

#     def update_transcripts(self, transcript):
#         if self.min_coord:
#             if transcript.start < self.min_coord:
#                 self.min_coord = transcript.start
#         else:
#             self.min_coord = transcript.start

#         if self.max_coord:
#             if transcript.end > self.max_coord:
#                 self.max_coord = transcript.end
#         else:
#             self.max_coord = transcript.end

#         if self.transcripts:
#             self.transcripts += [transcript]
#             self.chr += [transcript.chr]
#         else:
#             self.transcripts = [transcript]
#             self.chr = [transcript.chr]


####ENDS models
class Reader(object):
    def __init__(self, filename):
        self._filename = os.path.basename(filename)
        self._handler = open(filename, 'r')

    def __iter__(self):
        return self

    def next(self):
        try:
            line = next(self._handler)
        except ValueError:
            raise StopIteration
        return Record.fromline(line)

    # python 3 compatibility
    def __next__(self):
        return self.next()

    def close(self):
        self._handler.close()

class Record(object):
    def __init__(self, geneName, name, chrom, strand, txStart, txEnd,cdsStart, cdsEnd, exonCount, exonStarts, exonEnds):
        self._gene = geneName
        self._tr_name = name
        self._chr = chrom
        self._strand = strand
        self._tx_start = txStart
        self._tx_end = txEnd
        self._cds_start = cdsStart
        self._cds_end = cdsEnd
        self._exon_count = exonCount
        self._exon_start = exonStarts
        self._exon_ends = exonEnds

    @property
    def gene(self):
        return str(self._gene)

    @property
    def transcript(self):
        return str(self._tr_name)

    @property
    def chromosome(self):
        return str(self._chr)

    @property
    def strand(self):
        return str(self._strand)

    @property
    def txStart(self):
        return int(self._tx_start)

    @property
    def txEnd(self):
        return int(self._tx_end)

    @property
    def cdsStart(self):
        return int(self._cds_start)

    @property
    def cdsEnd(self):
        return int(self._cds_end)

    @property
    def n_exons(self):
        return int(self._exon_count)

    @property
    def exonStarts(self):
        return [int(x) for x in self._exon_start]

    @property
    def exonEnds(self):
        return [int(x) for x in self._exon_ends]

    @property
    def exons(self):
        return Exon.fromrecord(self)

    @property
    def cds_exons(self):
        """
        Return those exons which lie within the cds
        Also returns those partially inside the cds
        :return:
        """
        return [x for x in self.exons if x.stop >= self.cdsStart and x.start <= self.cdsEnd]

    def to_dict(self):
        d = {}
        d["geneName"] = self.gene
        d["name"] = self.transcript
        d["chrom"] = self.chromosome
        d["strand"] = self.strand
        d["txStart"] = self.txStart
        d["txEnd"] = self.txEnd
        d["cdsStart"] = self.cdsStart
        d["cdsEnd"] = self.cdsEnd
        d["exonStarts"] = self.exonStarts
        d["exonEnds"] = self.exonEnds
        d["exonCount"] = self.n_exons

        return d

    @property
    def line(self):
        line = []
        d = self.to_dict()
        for nlc in NUMERIC_LIST_COLUMNS:
            d[nlc] = ",".join(map(str, d[nlc])) + ","
        for col in COLUMNS:
            line += [d[col]]
        return "\t".join(map(str, line))



    @classmethod
    def fromdict(cls, items):
        """
        Builds a record from a dictionary.
        This dictionary must contain all fields specified in generics.COLUMNS
        """
        normal_columns = set(COLUMNS) - set(NUMERIC_LIST_COLUMNS) # <-- remember, this is UNORDERED!

        # first check whether all columns are there and properly formatted
        for c in COLUMNS:
            if c not in items:
                raise ValueError("Item {c} must be given".format(c=c))
        for nlc in NUMERIC_LIST_COLUMNS:
            if not isinstance(items[nlc], list):
                raise ValueError("Item {nlc} must be a list of integers".format(nlc=nlc))
            elif not all([isinstance(x, int) for x in items[nlc]]):
                raise ValueError("Item {nlc} must be a list of integers".format(nlc=nlc))
        for nc in NUMERIC_COLUMNS:
            if not isinstance(items[nc], int):
                raise ValueError("Item {nc} must be an integer".format(nc=nc))

        r = Record(**items)
        return r

    @classmethod
    def fromline(cls, line):
        """
        Builds a record from a line
        """
        raw_items = line.strip().split('\t')

        assert len(raw_items) >= 11, "Contains less than 11 columns!"
        items = dict()
        for i in range(11):
            items[COLUMNS[i]] = raw_items[i]
        for nc in NUMERIC_COLUMNS:
            items[nc] = int(items[nc])
        for lnc in NUMERIC_LIST_COLUMNS:
            if not items[lnc].endswith(','):
                raise ValueError("Malformed refFlat file! Value {lnc} must end in a comma".format(lnc=lnc))

            it = items[lnc].split(',')
            it.pop()
            items[lnc] = [int(x) for x in it]

        r = Record(**items)
        return r

if __name__ == "__main__": 

    reader = Reader(filename="refFlat.txt")
    print("exon")

    for record in reader:
        total_transcript_size=0
        for exon in record.exons:
            exon_size = exon.stop - exon.start #I dont substract 1 from this difference because the refflat file commonly is in 0 based coordinates at the start coordinates and 1-based at the end coordinates, Also I checked the coordinates in Genebank of multiple transcripts of the file and confirmed that this is the case for this file. 
            total_transcript_size += exon_size
        print("gene:%s total_transcript_size: %d" % ( exon.gene,total_transcript_size ))
