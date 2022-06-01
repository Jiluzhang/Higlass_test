## Pileups and average features
## https://cooltools.readthedocs.io/en/latest/notebooks/pileup_CTCF.html
## Pileups is beneficail for reliably observing features in low-coverage Hi-C or single-cell HiC maps
## On-diagonal pileup & Off-diagonal pileup

## import packages
import cooler
import bioframe
import numpy as np
import cooltools
import matplotlib.pyplot as plt

## Load data
clr = cooler.Cooler('W_rep1.multires.cool::/resolutions/10000')
resolution = clr.binsize

## Fetch and select genomic features
mm10_chromsizes = bioframe.fetch_chromsizes('mm10')
mm10_cens = bioframe.fetch_centromeres('mm10')  # empty DataFrame
mm10_arms = bioframe.make_chromarms(mm10_chromsizes, mm10_cens)
mm10_arms = mm10_arms[mm10_arms.chrom.isin(clr.chromnames)].reset_index(drop = True)
mm10_arms = mm10_arms[mm10_arms['chrom'] != 'chrM']  # chrM is too short

## Load features for anchors
h4k20me1 = bioframe.read_table('WH4K20_10kb.bedgraph')
h4k20me1.columns = ['chrom', 'start', 'end', 'score']
h4k20me1 = h4k20me1.query(f'chrom in {clr.chromnames}')
h4k20me1['mid'] = (h4k20me1.end + h4k20me1.start) // 2

## Select sites (top 20%)
sites = h4k20me1[h4k20me1['score'] >= np.percentile(h4k20me1['score'], 99)]

## On-diagonal pileup of observed over expected interactions
expected = cooltools.expected_cis(clr, view_df = mm10_arms, nproc = 20, chunksize = 1_000_000)
stack = cooltools.pileup(clr, sites, view_df = mm10_arms, expected_df = expected, flank = 300_000)
mtx = np.nanmean(stack, axis = 2)

## Plot heatmap
flank = 300_000
plt.imshow(
    np.log2(mtx),
    vmax = 0.25,
    vmin = -0.25,
    cmap='coolwarm',
    interpolation='none')

plt.colorbar(label = 'log2 mean obs/exp')
ticks_pixels = np.linspace(0, flank*2//resolution,5)
ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
plt.xticks(ticks_pixels, ticks_kbp)
plt.yticks(ticks_pixels, ticks_kbp)
plt.xlabel('relative position, kbp')
plt.ylabel('relative position, kbp')

plt.savefig('W_rep1_h4k20me1.pdf')
plt.close()


