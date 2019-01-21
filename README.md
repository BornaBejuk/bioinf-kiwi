# bioinf-kiwi

Installation and running:

git clone

g++ -std=c++0x *.cpp

./a.out pathCR data/EColi-synthetic/overlaps-c-r.paf pathRR data/EColi-synthetic/overlaps-r-r.paf pathFastaCtgs data/EColi-synthetic/ecoli_test_contigs.fasta pathFastaReads data/EColi-synthetic/ecoli_test_reads.fasta pathFastaOut data/EColi-synthetic/final.fasta SImin 0.9 maxDepth 40 nTimes 50

or 

./a.out running_example.txt
