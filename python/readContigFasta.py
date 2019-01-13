__author__ = 'Mario'

from proba_puta import *

def readFasta(filename, list_of_names):
    file = open(filename, "r")
    names = []
    contigs = []
    readNext=False
    for read in file:
        if read[0] == ">":
            name = read[1:].split()[0]
            if name in list_of_names:
                names.append(name.strip())
                readNext = True
        else:
            if readNext:
                contigs.append(read.strip())
                readNext = False
    file.close()
    return names, contigs


def readFastq(filename, list_of_names):
    file = open(filename, "r")

    reads = []
    names = []
    i = 0
    readNext = False
    for read in file:
        if i % 4 == 0:
            name = read[1:]
            if name in list_of_names:
                names.append(name.strip())
                readNext = True
        elif i % 4 == 1:
            if readNext:
                reads.append(read.strip())
                readNext = False
        i += 1

    file.close()
    return names, reads

def createReadingLists1(reads_filename, contigs_filename, list_of_names ):
    names, reads = [],[]
    n, r = readFasta(contigs_filename, list_of_names)
    names += n
    reads += r

    n, r = readFastq(reads_filename, list_of_names)
    names += n
    reads += r
    return names, reads

def createReadingLists(reads_filename, contigs_filename, list_of_names):
    names, reads = [],[]
    n, r = readFasta(contigs_filename, list_of_names)
    names += n
    reads += r

    n, r = readFasta(reads_filename, list_of_names)
    names += n
    reads += r
    return names, reads


def createFastaReference(connections, reads_filename, contigs_filename):

    results = find_connections(connections) #connections je bejukov ulaz tipa ['ctg1', 'read00001', 'read00007', 'read00002', 'ctg2']
    #print(results)
    names, reads = createReadingLists(reads_filename, contigs_filename, results)
    result_names = results # e sad njih treba postpajati. ali kako.
    result_reads = []

    #print(result_names)
    for i in result_names:
        result_reads.append(reads[names.index(i)]) ## e sad, njih treba spojiti, prema onim fjama iz utils-a.

    file = open("Reference.fasta", "w")
    title = ">Result\n"
    file.write(title)
    for i in result_reads:
        file.write(i) # to jest i.getContext() da se dobiju slova
    file.close()

def create_solution(reads_filename, contigs_filename, read_pairs):
    createFastaReference(read_pairs, reads_filename, contigs_filename)
    return


def example():
    reads = "C:\\Users\Mario\Desktop\Bioinformatika-scaffolding\EColi-synthetic\ecoli_test_reads.fasta"
    contigs = "C:\\Users\Mario\Desktop\Bioinformatika-scaffolding\EColi-synthetic\ecoli_test_contigs.fasta"
    input = [['ctg1', 'read00001', 'read00007', 'read00002', 'ctg2'], ["ctg3", 'read00002', 'read00004', 'read00006', "ctg1"], ["ctg3",  'read00001', 'read00007', 'read00002', "ctg2"]]
    create_solution(reads, contigs, input)

example()