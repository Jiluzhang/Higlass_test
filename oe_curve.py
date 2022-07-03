## Pileups and average features
## https://cooltools.readthedocs.io/en/latest/notebooks/pileup_CTCF.html
## Pileups is beneficail for reliably observing features in low-coverage Hi-C or single-cell HiC maps
## On-diagonal pileup & Off-diagonal pileup

## import packages
import cooler
import bioframe
import pandas as pd
import numpy as np
import cooltools
import matplotlib.pyplot as plt
import argparse
import warnings
from plotnine import *

warnings.filterwarnings('ignore')

#bedtools makewindows -g ~/ref_genome/chrom_sizes/mm10.chrom.sizes -w 10000 > mm10_10kb.bed
#bedtools intersect -a mm10_10kb.bed -b WH4K20_10kb.bedgraph -wa -wb | awk '{print $1 "\t" $2 "\t" $3 "\t" $7}' | \
#awk '{if($3%10000==0) print$0}' > WH4K20_10kb_all.bedgraph   # pay attention to end of chromosome

resolution = 10_000

## genomic features
mm10_arms = pd.read_csv('mm10_arms.csv')
chr = ['chr' + str(i) for i in range(1, 20)]
chr.append('chrX')

## Load features for anchors
h4k20me1 = bioframe.read_table('WH4K20_10kb_all.bedgraph')
h4k20me1.columns = ['chrom', 'start', 'end', 'score']
h4k20me1 = h4k20me1.query(f'chrom in {chr}')
h4k20me1['mid'] = (h4k20me1.end + h4k20me1.start) // 2
h4k20me1 = h4k20me1[h4k20me1['score'] !=0]

## Select sites (high & low)
hig = h4k20me1[h4k20me1['score'] >= np.percentile(h4k20me1['score'], 70)]
low = h4k20me1[h4k20me1['score'] <  np.percentile(h4k20me1['score'], 70)]
# 0.367184

dict = {'high': hig, 'low': low}

def plot_oe_vs_dist(c, t):
    wt_clr = cooler.Cooler(c + '.multires.cool::/resolutions/10000')
    wt_expected = cooltools.expected_cis(wt_clr, view_df = mm10_arms, nproc = 40, chunksize = 1_000_000)
    
    mt_clr = cooler.Cooler(t + '.multires.cool::/resolutions/10000')
    mt_expected = cooltools.expected_cis(mt_clr, view_df = mm10_arms, nproc = 40, chunksize = 1_000_000)
    
    for lev in ['high', 'low']:
        paired_sites = bioframe.pair_by_distance(dict[lev], min_sep = 20_000, max_sep = 1_000_000, suffixes = ('1', '2'))
        paired_sites.loc[:, 'mid1'] = (paired_sites['start1'] + paired_sites['end1'])//2
        paired_sites.loc[:, 'mid2'] = (paired_sites['start2'] + paired_sites['end2'])//2
        paired_sites['dist'] = abs(paired_sites['mid1'] - paired_sites['mid2'])
        
        wt_stack = cooltools.pileup(wt_clr, paired_sites, view_df = mm10_arms, expected_df = wt_expected, flank = 0, nproc = 40)
        wt_tmp = paired_sites[['dist']]
        wt_tmp['oe'] = wt_stack[0][0]
        wt_res = wt_tmp[['dist', 'oe']].groupby('dist', as_index = False).mean()
        wt_res['idx'] = 'wt'
        
        mt_stack = cooltools.pileup(mt_clr, paired_sites, view_df = mm10_arms, expected_df = mt_expected, flank = 0, nproc = 40)
        mt_tmp = paired_sites[['dist']]
        mt_tmp['oe'] = mt_stack[0][0]
        mt_res = mt_tmp[['dist', 'oe']].groupby('dist', as_index = False).mean()
        mt_res['idx'] = 'mt'
        
        res = pd.concat([wt_res, mt_res])
        res['dist'] = res['dist'] / 1000000
        
        ## Plot heatmap
        p = ggplot(res, aes(x = 'dist', y = 'oe', color = 'idx')) + geom_smooth() + \
            scale_x_continuous(limits = [0, 1], breaks = np.arange(0, 1+0.01, 0.1)) + \
            scale_y_continuous(limits = [0, 1.5], breaks = np.arange(0, 1.5+0.3, 0.3)) + theme_bw()
        p.save(c + '_' + t + '_h4k20me1_' + lev + '_oe_vs_dist.pdf')
        
        print(c, t, lev, 'h4k20me1 done')

plot_oe_vs_dist(c = 'W1G1', t = 'R1G1')
