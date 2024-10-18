# Mission 3

1. For each of 22 autosomes and 2 sex chromosomes in human, read an hg38 unmasked sequence and store it into a single string.
2. Read the RefFlat file, and parse the information, which you have done in Mission 2 already.
3. For each RefSeq sequence, read the genomic sequence of each exon, partition them into 5’UTR, CDS, and 3’UTR by combining 1 and 2.

```diff
[Data Structure]
    RefSeq ID
    Gene Symbol
    Chromosome ID	(e.g., 1,2,…,22,23 for X, and 24 for Y)
    Strand
    Num of Exons
    Exon Start Positions,		Exon End Positions
-   Exon Seq Size			Exonic Genome Sequence
-   5’UTR Seq Size			5’UTR Genome Sequence
-   ORF Seq Size			ORF Genome Sequence
-   3’UTR Seq Size			3’UTR Genome Sequence
```

## More Details on Mission 3

1. Use the RefFlat file provided on 08/05/18. Count the number of entries in the RefFlat file.
2. Remove non-NM sequences and keep the ones aligned only on chr 1-22, X, or Y. Count the number of entries left.
3. Remove NM sequences that have multiple entries in the RefFlat file. Count the number of entries left.
4. Remove NM sequences that have wrong ORFs: no start codon, no stop codon, ORF size of 3N+1 or 3N+2, or internal stop codons. Count the number of entries left.
5. For each gene that has multiple isoforms, select a representative isoform by picking up a RefSeq sequence that has the lowest NM ID(an integer). Count the number of entries left.
6. Sort the list of the RefSeq class objects by the ascending order of the RefSeq ID (integer), and store the RefSeq ID, the gene symbol, 5’UTR size, ORF size, and 3’UTR size in the provided template file. When computing ORF size, include the terminal stop codon (3nt).
7. When storing ORF sequences, include both start and stop codons.

Estimated run time: < 3 mins
