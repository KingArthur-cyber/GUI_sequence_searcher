Input
It takes a fasta file as an input and searches those sequences within a biggeer genome.
Output:
Output is a text file pointinf to the location of the subsequence within the sequence.
Working:
It is based upon BLAST's method of forming kmers, and then aligning them with the refseq, you could use various algorithms to judge the alignment, like edit distance, needleman wunsch algorithm, the option to choose algorithm for alignment 
is underway...
