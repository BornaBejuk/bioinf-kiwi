__author__ = 'Mario'

from proba_puta import *

def readFasta(filename):
    file = open(filename, "r")
    names = []
    contigs = []
    for read in file:
        if read[0] == ">":
            name = read[1:].split()[0]
            names.append(name.strip())
        else:
            contigs.append(read.strip())
    file.close()
    return names, contigs


def readFastq(filename):
    file = open(filename, "r")

    reads = []
    names = []
    i = 0
    for read in file:
        if i % 4 == 0:
            name = read[1:]
            names.append(name.strip())
        elif i % 4 == 1:
            reads.append(read.strip())
        i += 1

    file.close()
    return names, reads

def createReadingLists1(reads_filename, contigs_filename):
    names, reads = [],[]
    n, r = readFasta(contigs_filename)
    names += n
    reads += r

    n, r = readFastq(reads_filename)
    names += n
    reads += r
    return names, reads

def createReadingLists(reads_filename, contigs_filename):
    names, reads = [],[]
    n, r = readFasta(contigs_filename)
    names += n
    reads += r

    n, r = readFasta(reads_filename)
    names += n
    reads += r
    return names, reads


def createFastaReference(read_pairs, names, reads):

    results = connect_into_fasta(read_pairs)

    result_names = results[0] # ili neki drugi odabir
    result_reads = []
    print(result_names)
    print(names)
    for i in result_names:
        result_reads.append(reads[names.index(i)])

    file = open("Reference.fasta", "w")
    title = ">Result\n"
    file.write(title)
    for i in result_reads:
        file.write(i) # to jest i.getContext() da se dobiju slova
    file.close()

def create_solution(reads_filename, contigs_filename, read_pairs):
    names, reads = createReadingLists(reads_filename, contigs_filename)
    createFastaReference(read_pairs, names, reads)
    return


def example():
    reads = "C:\\Users\Mario\Desktop\Bioinformatika-scaffolding\EColi-synthetic\ecoli_test_reads.fasta"
    contigs = "C:\\Users\Mario\Desktop\Bioinformatika-scaffolding\EColi-synthetic\ecoli_test_contigs.fasta"
    pairs = [("read00001", "ctg1"), ("ctg2", "ctg1"), ("ctg3", "ctg2"), ("ctg3", "ctg1")]
    create_solution(reads, contigs, pairs)

example()