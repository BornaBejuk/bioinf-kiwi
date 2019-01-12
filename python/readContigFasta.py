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


def createFastaReference(connections, names, reads):

    results = find_connections(connections) #connections je bejukov ulaz tipa ['ctg1', 'read00001', 'read00007', 'read00002', 'ctg2']
    #print(results)
    result_names = results # ili neki drugi odabir
    result_reads = []
    #print(result_names)
    #print(names)
    for i in result_names:
        result_reads.append(reads[names.index(i)]) ## e sad, njih treba spojiti, prema onim fjama iz utils-a.

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
    input = [['ctg1', 'read00001', 'read00007', 'read00002', 'ctg2'],
             ["ctg3", 'read00002', 'read00004', 'read00006', "ctg1"], ["ctg3",  'read00001', 'read00007', 'read00002', "ctg2"]]
    create_solution(reads, contigs, input)

example()