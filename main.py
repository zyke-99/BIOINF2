from random import sample

from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


class FastqRecord:
    def __init__(self, id, sequence, third, quality):
        self.id = id
        self.sequence = sequence
        self.third = third
        self.quality = quality


class OccurRecord:
    def __init__(self, id, sequence, gc):
        self.id = id
        self.sequence = sequence
        self.gc = gc


enc = {
    'Sanger Phred+33': [33, 73],
    'Solexa Solexa+64': [59, 104],
    'Illumina 1.3+ Phred+64': [64, 126],
    'Illumina 1.5+ Phred+64': [67, 105],
    'Illumina 1.8+ Phred+33': [33, 74]
}


def determine_encoding(unique_characters):
    potential_encodings = []
    for encoding in enc:
        is_of_encoding = True
        for character in unique_characters:
            if ord(character) < enc[encoding][0] or ord(character) > enc[encoding][1]:
                is_of_encoding = False
        if is_of_encoding:
            potential_encodings.append(encoding)
    return potential_encodings


def read():
    return SeqIO.parse("reads_for_analysis.fastq", "fastq")


def read2():
    records = []
    with open("reads_for_analysis.fastq", "r") as f:
        id = f.readline().strip()
        while id:
            if id:
                sequence = f.readline().strip()
                third = f.readline().strip()
                quality = f.readline().strip()
                records.append(FastqRecord(id, sequence, third, quality))
                id = f.readline().strip()
        return records


def get_all_unique_chars(string):
    return set(string)


def plot(occur_records):
    xpoints = []
    ypoints = []
    for x in range(0, 100):
        xpoints.append(x/100)
        ypoints.append(len(list(filter(lambda record: record.gc == x/100, occur_records))))
    plt.plot(xpoints, ypoints, 'o')
    plt.show()


def blast_search(seq):
    result_handle = NCBIWWW.qblast("blastn", "nt", seq, alignments=1, hitlist_size=1)
    records = NCBIXML.parse(result_handle)
    for rec in records:
        for alignment in rec.alignments:
            return alignment.hit_def


def calc_occurs(records):
    occurs = []
    for record in records:
        gc = (record.sequence.count('G') + record.sequence.count('C')) / len(record.sequence)
        gc = round(gc, 2)
        occurs.append(OccurRecord(record.id, record.sequence, gc))
    return occurs


def get_records_from_peaks(occur_records):
    records = []
    records += (list(filter(lambda record: record.gc == 0.54, occur_records)))[0:5]
    records += (list(filter(lambda record: record.gc == 0.34, occur_records)))[0:5]
    records += (list(filter(lambda record: record.gc == 0.70, occur_records)))[0:5]
    return records


if __name__ == '__main__':
    conc_str = ''
    records = read2()
    for record in records:
        conc_str += record.quality
    print(determine_encoding(get_all_unique_chars(conc_str)))
    occurs = calc_occurs(records)
    plot(occurs)
    bacterias = []
    for record in get_records_from_peaks(occurs):
        bacteria = blast_search(record.sequence)
        print(record.id)
        print(bacteria)
        bacterias.append(bacteria)
