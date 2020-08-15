# Naive_exact_matching
Context: By using and modifying Python programs provided in the Coursera course entitled 'Algorithms for DNA Sequencing', I checked the quality of reads from a DNA sequencing experiment and aligned reads with a genome.

What does the program do?
1. readGenome turns the content of a FASTA file into a string. This was provided for the assignment.
2. readFastq turns the content of a FASTQ file into two parallel lists containing the reads and the base quality scores of each read respectively. This was provided for the assignment.
3. reverseComplement turns a string representing a DNA sequence into its reverse complement. This was provided for the assignment.
4. naive returns the offset of each alignment where read p matches genome t. It is the naive exact matching algorithm. This was provided for the assignment.
5. naive_with_rc is a modified version of naive. I coded it to also consider the reverse complement of read p in the alignment task.
6. naive_2mm is another modified version of naive. I coded it to accept alignments with up to two mismatches. Therefore, it is no longer an exact matching algorithm.
7. phred33ToQ converts an ASCII character to a base quality score. I coded it.
8. The final block of code extract the reads and the base quality scores of each read from a FASTQ file, calculates the number of reads and the length of each read, calculates the average quality score at each offset (in each sequencing cycle), and identifies the offset with the lowest average quality score. I coded it.

Software: Python 3.7.

Input files:
1. A FASTA file containing a genome. The lambda virus genome is included as an example (lambda_virus.fa).
2. A FASTQ file. A file containing real DNA sequencing reads derived from a human is included as an example (ERR037900_1.first1000.fastq).
