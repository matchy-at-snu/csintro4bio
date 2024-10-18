# Mission 5: Predicting the Overexpressed miRNAs

Repeat the analysis given in Mission 4 for each the provided 3 transcriptome data. The analyses are identical to those in Mission 4 except for the following changes:

Correct your P values by using the bonferonni correction. When counting the number of tests for the bonferonni correction, only consider Fisher’s exact tests whose relative risk are bigger than 1.0 (A/B>1.0).
Report the 5 most significantly enriched motifs in the highly downregulated genes, in the ascending order of P values.
For those highly associated motifs, report whether they are 7mer-m8 or 7mer-A1 sites of known human miRNAs.
Report your results in the provided Mission5_template.

Repeat the above analyses for ORF 7mers, and discuss about your results in comparison with 3’UTR results and about possible methodological improvements.

Estimated run time: <60 mins

## Result discussion

"According to the analysis of the 3'UTR and ORF miRNA targetting sites, it is easy to see that the top 5 motifs for 3'UTR had a very low p-value (< 0.01) even after the bonferonni correction. In the case of ORFs, the results were bad: the analysis generates no significant results. For datasset1, even though the top 5 motifs all had a p value lower than 0.01, none of them matched with know miRNAs. In dataset2 and dataset 3, the adjusted p-values were so large, with some even reaching 1. One noteworthy fact is that dataset2 has the same top1 motif in both ORF and 3'UTR group, but the p-value in ORF group is higher to the extent of losing significance and the relative risk (degree of enrichment) is lower. Since bonferroni method is very stringent, alternative analyses using Benjamini-Hochberg procedure to control false discovery rate were performed. However, the results were similar: the analysis on the motifs in ORF region no significant (p-value cutline 0.01) results. However, if we set the cutline to 0.05, there is one significant generated in dataset 2. In addition, the analyses of the above procedure were also performed for 5'UTR. As a result, almost all the motifs had very high p-values, and those few with lower p-values have no matched miRNAs, thus suggesting no significant results.

From the results it is easy to see that the miRNA rugulatory efficiency is relatively low in 5'UTR and ORF region compared to 3'UTR. It's very likely that the steric hinderance of the 5'cap and the translation machinery make silencing less efficient in 5'UTR and ORF region. Also, recent research showed that the interaction between binding sites in coding region and 3'UTRs might make the effect of miRNA regulation milder, therefore miRNA target sites in coding regions are under negative selection[1]. This might also explain the scarcity of siginificant results in ORF region.

In the above analysis, when constructing the contingency table for the fisher's exact test, the fold-change cutoff value of highly downregulated was set to -0.5. This criterion might be a little bit too loose since many motifs had very low there seemed to be multiple significant results. It might be okay to set a tighter cirterion like -1 or lower.

1.Fang, Z., & Rajewsky, N. (2011). The impact of miRNA target sites in coding sequences and in 3′ UTRs. PloS one, 6(3), e18067."
