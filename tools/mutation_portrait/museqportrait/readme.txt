This is a package to generate post-processing reports ('portraits') from the mutationSeq vcf output.

v0.99.0 (released December 19, 2013):

Plot types:

1. Probability density plots

2. Bean plots of allele ratios across a range of P-values:
- mutations are binned by P-value such that a plot labeled "0.5" on the x-axis contains mutations with P-values in the range of 0.5-0.6, and so on. The estimated density of mutations across allele ratios is plotted as a gray curve. Individual mutation values are shown as horizontal black lines. The thick and wider black line indicates the mean.

3. Coverage versus allele ratio plots:
- the estimated density is plotted 
- individual points are too numerous to show
- x and y axes are scaled to capture +/- 3 standard deviations

4. Tri-nucleotide context distribution
- Inspired by plots in Alexandrov et al. Cell Rep. 2013 Jan 31;3(1):246-59 (PMID: 23318258)

5. Genome-wide mutation distribution
- Genome length is divided into 3000 bins. First, each chromosome is assigned a number of bins proportional to its size (i.e. largest chromosome assigned the most bins, etc.). Second, the bin size is computed for each chromosome as the total chromosome size divided by this bin number. Third, mutations are enumerated for each bin and then divided by the bin size in Mbp. This approach ensures that the bin sizes are the same for any one chromosome.
