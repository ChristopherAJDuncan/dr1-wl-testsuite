# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2025-06-24 08:54:57
# @Last Modified by:   lshuns
# @Last Modified time: 2025-12-28 18:35:08

### With the given split criteria, compare key properties of subsamples
###### red vs blue galaxies

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

import plotting

## >>>>>>>>>>>>> I/O

## which shear measurement
shear_type = 'lensmc'
# shear_type = 'metacal'

## u-r cut value (selection criteria)
cut_value = 1.2

## number of tomo bins
Ntomo = 6

## where to find catalogues
TITLE = 'euclid_tr1_a_v1 + euclid_tr1_b_v1_1'
inpath_list = ['/sdf/data/kipac/u/liss/euclid/TR1/23538.parquet', 
               '/sdf/data/kipac/u/liss/euclid/TR1/23539.parquet']
used_cols_list = [['object_id', 
                   'right_ascension', 
                   'declination', 
                   'flux_radius',
                   'phz_pp_mode_mu', 
                   'phz_pp_mode_mr', 
                   'sersic_sersic_vis_radius', 
                   'sersic_sersic_vis_index', 
                   'sersic_sersic_vis_axis_ratio',
                   'flux_detection_total'],
                  ['object_id', 
                   'tom_bin_id', 
                   'phz_mode_1', 
                   f'she_{shear_type}_weight']]

## which columns to plot
plot_cols = ['flux_radius', 
             'phz_mode_1', 
             'sersic_sersic_vis_radius', 
             'sersic_sersic_vis_index', 
             'sersic_sersic_vis_axis_ratio',
            'mag_detection_total']
XRANGE_list = [[1, 20], 
               [0, 2.5],  
               [0.01, 20],
               [0, 9],
               [0, 1],
               [20, 25]]
xlog_list = [True, 
             False,
             True,
             False,
             False,
             False
             ]

## >>>>>>>>>>>>>> Workhorse

## +++ Load and merge catalogue
for i_cata, inpath in enumerate(inpath_list):
    if i_cata == 0:
        cata = pd.read_parquet(inpath)[used_cols_list[i_cata]]
        print(f">>> Number loaded from cata {i_cata}", len(cata))
    else:
        cata_tmp = pd.read_parquet(inpath)[used_cols_list[i_cata]]
        print(f">>> Number loaded from cata {i_cata}", len(cata_tmp))
        cata = cata.merge(cata_tmp, on='object_id')
        del cata_tmp
        print(f">>>>>> Number after merging", len(cata))
## weight selection
cata = cata.loc[(cata[f'she_{shear_type}_weight']>0)].reset_index(drop=True)
print(f">>> Number after she_{shear_type}_weight selection", len(cata))
## additional data quality filters
cata = cata[(cata['flux_detection_total']>0)].reset_index(drop=True)
print(f">>> Number after flux quality filter", len(cata))
## calculate magnitude
cata['mag_detection_total'] = 23.9 - 2.5 * np.log10(cata['flux_detection_total'].values)

## +++ Assign subsample ID
cata['subsample_ID'] = 0
## calculate colour and cut
cata.loc[((cata['phz_pp_mode_mu'].values - cata['phz_pp_mode_mr'].values)>=cut_value), 'subsample_ID'] = 1
print(f">>> subsamples defined with cut_value = {cut_value}")
print(">>>>>>> subsample number ratio (blue/red)", np.sum(cata['subsample_ID']==0)/np.sum(cata['subsample_ID']==1))
print(">>>>>>> subsample weight ratio (blue/red)", np.sum(cata.loc[cata['subsample_ID']==0, f'she_{shear_type}_weight'].values)/np.sum(cata.loc[cata['subsample_ID']==1, f'she_{shear_type}_weight'].values))

## +++ 1D hist without tomo bins
mask_red = (cata['subsample_ID']==1)
for i_col, plot_col in enumerate(plot_cols):
    outpath = f'./plots/B_red_blue_1Dhist_all_{plot_col}_{shear_type}_weighted.png'
    XLABEL = plot_col
    xlog = xlog_list[i_col]
    XRANGE = XRANGE_list[i_col]
    paras = [cata.loc[mask_red, plot_col].values, 
                cata.loc[(~mask_red), plot_col].values]
    wgs = [cata.loc[mask_red, f'she_{shear_type}_weight'].values, 
                cata.loc[(~mask_red), f'she_{shear_type}_weight'].values]
    COLORs = ['r', 'b']
    LABELs = [f'u-r>={cut_value}', f'u-r<{cut_value}']
    nbins = 60
    DENSITY = False
    YLABEL = f'{shear_type}-weighted counts'
    plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=TITLE, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, ylog=False,
                loc_legend='best', 
                LABEL_position='inSub', LABEL_cols=1,
                font_size=12, usetex=False,
                cumulative=False, 
                FIGSIZE=[6.4, 4.8],
                LINEs=None, LINEWs=None,
                TIGHT=False,
                alpha=None,
                HISTTYPE_list=None)

## +++ 1D hist with tomo bins
N_plots = Ntomo
mask_red = (cata['subsample_ID']==1)
for i_col, plot_col in enumerate(plot_cols):
    outpath = f'./plots/B_red_blue_1Dhist_tomo_{plot_col}_{shear_type}_weighted.png'
    XLABEL = plot_col
    xlog = xlog_list[i_col]
    XRANGE = XRANGE_list[i_col]
    paras_list = []
    wgs_list = []
    subLABEL_list = []
    for i_ZBbin in range(Ntomo):
        mask_tomo = cata['tom_bin_id']==(i_ZBbin+1)
        paras_list.append([cata.loc[mask_red&mask_tomo, plot_col].values, 
                           cata.loc[(~mask_red)&mask_tomo, plot_col].values])
        wgs_list.append([cata.loc[mask_red&mask_tomo, f'she_{shear_type}_weight'].values, 
                           cata.loc[(~mask_red)&mask_tomo, f'she_{shear_type}_weight'].values])
        subLABEL_list.append(f'tomo {i_ZBbin+1}')
    COLORs_list = [['r', 'b']] * N_plots
    LABELs_list = [[f'u-r>={cut_value}', f'u-r<{cut_value}']] * N_plots
    nbins_list = [60] * N_plots
    DENSITY = False
    YLABEL = f'{shear_type}-weighted counts'
    plotting.HistPlotFunc_subplots(outpath, N_plots,
                                paras_list, wgs_list, COLORs_list, LABELs_list,
                                nbins_list, XRANGE, YRANGE=None,
                                LINEs_list=None, LINEWs_list=None,
                                subLABEL_list=subLABEL_list, subLABEL_locX=0.1, subLABEL_locY=0.8,
                                XLABEL=XLABEL, YLABEL=YLABEL,
                                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                                xtick_spe=None, ytick_spe=None,
                                vlines_list=None, 
                                vline_styles_list=None, vline_colors_list=None, vline_labels_list=None, vline_widths_list=None,
                                hlines_list=None, 
                                hline_styles_list=None, hline_colors_list=None, hline_labels_list=None, hline_widths_list=None,
                                xlog=xlog, ylog=False,
                                loc_legend='best',
                                font_size=12, usetex=False,
                                LABEL_position='top', LABEL_position_SUBid=0,
                                LABEL_cols=2, 
                                FIGSIZE=[6.4, 4.8],
                                TIGHT=False,
                                HISTTYPEs_list=None,
                                shareX=True,
                                shareY=True,
                                XRANGE_list=None,
                                YRANGE_list=None)

## +++ Sky distribution
N_sky_bins = 60
ra_edges = np.unique(
    np.quantile(cata['right_ascension'], 
                q=np.linspace(0, 1, N_sky_bins + 1)))
dec_edges = np.unique(
    np.quantile(cata['declination'], 
                q=np.linspace(0, 1, N_sky_bins + 1)))
cata['ra_bin'] = pd.cut(cata['right_ascension'], 
                        bins=ra_edges, 
                        labels=False, 
                        include_lowest=True)
cata['dec_bin'] = pd.cut(cata['declination'], 
                         bins=dec_edges, 
                         labels=False, 
                         include_lowest=True)
## compute weighted mean RA and Dec per bin
def compute_weighted_center(g):
    weights = g[f'she_{shear_type}_weight']
    if len(g) == 0 or weights.sum() == 0:
        return pd.Series({'ra_center': np.nan, 'dec_center': np.nan})
    return pd.Series({
        'ra_center': np.average(g['right_ascension'], weights=weights),
        'dec_center': np.average(g['declination'], weights=weights),
    })
weighted_centers = cata.groupby(['ra_bin', 
                                 'dec_bin']).apply(
                                     compute_weighted_center,
                                     include_groups=False).reset_index()
weighted_centers = weighted_centers.dropna() 
## group by bin and subsample, sum weights
weight_sum = cata.groupby(['ra_bin', 
                           'dec_bin', 
                           'subsample_ID'])[f'she_{shear_type}_weight'].sum().unstack(fill_value=0)
weight_sum.columns = ['w_0', 'w_1']
## add small epsilon to avoid division by zero
weight_sum['ratio'] = weight_sum['w_1'] / (weight_sum['w_0'] + 1e-10)
## merge with sky info
plot_df = pd.merge(weight_sum.reset_index(), 
                    weighted_centers, 
                    on=['ra_bin', 'dec_bin'])
plot_df = plot_df.dropna()
## plot
plt.figure(figsize=(8, 6))
sc = plt.scatter(plot_df['ra_center'], 
                 plot_df['dec_center'], 
                 c=plot_df['ratio'],
                cmap='coolwarm', s=2,
                vmin=0.8, vmax=1.2)
plt.gca().invert_xaxis()
plt.xlabel('Weighted mean RA')
plt.ylabel('Weighted mean dec')
plt.title(TITLE)
plt.colorbar(sc, label=f'{shear_type}-weighted ratio (red / blue)')
plt.tight_layout()
plt.savefig(f'./plots/B_red_blue_sky_{shear_type}_weighted.png', dpi=300)