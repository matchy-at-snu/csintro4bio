# Mission 4: Finding Overrepresented Motifs

Using the provided log2(fold-change) values, discover the motifs that are significantly enriched in the highly downregulated genes after miR-1 transfection into HeLa cell.

Use:

1. the Fisher’s exact test for the $2 \times 2$ contingency table shown below,
2. a fold-change cutoff of $-0.50$ to decide the highly downregulated genes (not including $-0.50$), and
3. a motif size of $7$, which means that we will iterate the test `4**7=16,384` times and pick the most significantly enriched motifs.

|                 | Motif | NotMotif |
|-----------------|---------------|-------------------|
| HighlyDownregulated | $n_1$             | $n_2$                 |
| NotHighlyDownregulated  | $n_3$             | $n_4$                 |

## Relative Risk v.s. Odds Ratio

|                 | Disease | Not Disease |
|-----------------|---------------|-------------------|
| Exposure | $n_1$             | $n_2$                 |
| Unexposure  | $n_3$             | $n_4$                 |

### Relative task 
$$
A = n_1 / (n_1+n_2)
B = n_3 / (n_3+n_4)
$$

$A / B$: relative risk.

### Disease odds ratio

- Odds in favor of disease in exposed group:		$C = n_1 / n_2$
- Odds in favor of disease in unexposed group:	$D = n_3 / n_4$

$$
C / D = \frac{n_1 / n_2}{n_3 / n_4} = \frac{n_1 n_4}{n_2 n_3}
$$

- Odds ratio > 1: a greater likelihood of disease among the exposed.
- Odds ratio < 1: a greater likelihood of disease among the unexposed.

## Additional analysis

Compute the fraction of genes that include the motif in highly downregulated genes (A) and the fraction of genes that include the motif in those that are not highly downregulated (B).

Then, the enrichment of the motif is $A/B$, and we will pick the motifs that have $A/B > 1.0$.

Report the 10 most significantly enriched motifs in the highly downregulated genes, in the ascending order of P values using the provided template file.
