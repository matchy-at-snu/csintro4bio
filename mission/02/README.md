# Mission 2

1. Read the RefFlat file, and parse the information, and store them into a list of class objects with your own design of instance variables and methods.

2. Your class should include instance variables that correspond to the following data structure.

```toml
[Data Structure]
RefSeq ID
Gene Symbol
Chromosome ID (e.g., 1,2,…,22,23 for X, and 24 for Y)
Strand
Num of Exons
Exon Start Positions,
Exon End Positions
```

## More Details on Mission 2

1. Use the RefFlat file provided on 11/25/18. Count the number of entries in the RefFlat file.

2. Remove non-NM sequences and keep the ones aligned only on chr 1-22, X, or Y. Count the number of entries left.
3. Remove NM sequences that have multiple entries in the RefFlat file. Count the number of entries left.
4. Sort the list of the RefSeq class objects by the ascending order of the RefSeq ID (not a string but an integer), and store the RefSeq ID and the Gene Symbol in the provided template file.

## Writeup

The CG suppression seems to be common in mammalian genomes as well as bird genomes, and it is usually explained through the methylation-deamination-mutation hypothesis, whereby the methylation of CG to 5-methylcytosine (forming the famous CpG islands) and subsequent deamination to thymine results, if unrepaired, in conversion of CG to TG/CA. The methylation hypothesis is supported by the fact that invertebrates that do not possess the canonical methylases (e.g. Drosophila and C. elegans) do not exhibit significant CG dinucleotide suppression. Also, the TG/CA contents are also slightly increased in chicken and human.[1]

TA is the least stable dinucleotide stacking pair and is prominent in some regulatory signals (TATA box, 3' polyadenylation...). Avoidance of spurious signal sequences and considerations of DNA stability could both act to suppresss overall levels of TA. It is also possible that the TA dinucleotide is unfavored in coding region since UA is relatively susceptible to cleavage by ribonucleases in mRNA, but this preposition requires further examination.[1]

Other discrepancies between relative frequency and expected frequency mostly have something to do with nucleosome and isochore organization. According to recent research, there is a strong periodicity of AA and TT dinucleotides in C. elegans with a interval of ~10 nt. It is believed to demarcate the nucleosome regions[2]. Since the chromosomes of nearly all eukaryotes are packaged into nucleosomes, the flanking sequences should naturally be enriched. The dinucleotides flanking nucleosomes in human [CC/GG](3) are also slightly enriched. The level of the dinucleotides that belong to the isochore region, for example AT/CG periodicity in chicken[4], are also slightly elevated in their genome. In D. melanogaster, there is a large-scale periodic enrichment of AA/TT associated with nucleosome occupancy, and GC dinucleotide peaks in linker regions[5], which corresponds with the elevation in frequency of these dinucleotides.
