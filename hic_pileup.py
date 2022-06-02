## Pileups and average features
## https://cooltools.readthedocs.io/en/latest/notebooks/pileup_CTCF.html
## Pileups is beneficail for reliably observing features in low-coverage Hi-C or single-cell HiC maps
## On-diagonal pileup & Off-diagonal pileup

############# h4k20me1_hic_pileup.py ################
## import packages
import cooler
import bioframe
import pandas as pd
import numpy as np
import cooltools
import matplotlib.pyplot as plt
import argparse
import warnings

warnings.filterwarnings('ignore')

## parse parameters
parser = argparse.ArgumentParser(description = 'Hi-C pileup heatmap')
parser.add_argument('-S', '--sample', help = 'sample name')
args = parser.parse_args()
sample = args.sample

## Load data
clr = cooler.Cooler(sample + '.multires.cool::/resolutions/10000')
resolution = 10_000

## Fetch and select genomic features (need internet)
#mm10_chromsizes = bioframe.fetch_chromsizes('mm10')
#mm10_cens = bioframe.fetch_centromeres('mm10')  # empty DataFrame
#mm10_arms = bioframe.make_chromarms(mm10_chromsizes, mm10_cens)
#mm10_arms = mm10_arms[mm10_arms.chrom.isin(clr.chromnames)].reset_index(drop = True)
#mm10_arms = mm10_arms[mm10_arms['chrom'] != 'chrM']  # chrM is too short
#mm10_arms.to_csv('mm10_arms.csv', index = False)
mm10_arms = pd.read_csv('mm10_arms.csv')
#chrom,start,end,name
#chr1,0,195471971,chr1_p
#chr2,0,182113224,chr2_p
#chr3,0,160039680,chr3_p
#chr4,0,156508116,chr4_p
#chr5,0,151834684,chr5_p


## Load features for anchors
h4k20me1 = bioframe.read_table('WH4K20_10kb.bedgraph')
h4k20me1.columns = ['chrom', 'start', 'end', 'score']
h4k20me1 = h4k20me1.query(f'chrom in {clr.chromnames}')
h4k20me1['mid'] = (h4k20me1.end + h4k20me1.start) // 2

## Select sites (high & mid & low)
hig = h4k20me1[h4k20me1['score'] >= np.percentile(h4k20me1['score'], 95)]
mid = h4k20me1[(h4k20me1['score'] >= np.percentile(h4k20me1['score'], 47.5)) &
               (h4k20me1['score'] <= np.percentile(h4k20me1['score'], 52.5))]
low = h4k20me1[h4k20me1['score'] <= np.percentile(h4k20me1['score'], 5)]
## hig: 1.29~10.34  mid: 0.21~0.24  low: 0~0.04
#sites = bioframe.cluster(sites, min_dist = resolution * 5).drop_duplicates('cluster').reset_index(drop = True)
# results from filteration and unfilteration are almost identity


## On-diagonal pileup of observed over expected interactions
expected = cooltools.expected_cis(clr, view_df = mm10_arms, nproc = 20, chunksize = 1_000_000)

dict = {'high': hig, 'middle': mid, 'low': low}

for lev in ['high', 'middle', 'low']:
    stack = cooltools.pileup(clr, dict[lev], view_df = mm10_arms, expected_df = expected, flank = 500_000, nproc = 20)
    mtx = np.nanmean(stack, axis = 2)
    
    ## Plot heatmap
    flank = 500_000
    plt.imshow(np.log2(mtx), vmax = -0.1, vmin = 0.1, cmap = 'coolwarm', interpolation='none')
    plt.colorbar(label = 'log2 mean obs/exp')
    ticks_pixels = np.linspace(0, flank*2//resolution,5)
    ticks_kbp = ((ticks_pixels-ticks_pixels[-1]/2)*resolution//1000).astype(int)
    plt.xticks(ticks_pixels, ticks_kbp)
    plt.yticks(ticks_pixels, ticks_kbp)
    plt.xlabel('relative position, kbp')
    plt.ylabel('relative position, kbp')
    plt.savefig(sample + '_h4k20me1_' + lev + '_pileup.pdf')
    plt.close()
    
    print(sample, lev, 'done')
