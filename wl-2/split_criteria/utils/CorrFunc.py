# -*- coding: utf-8 -*-
# @Author: ssli
# @Date:   2022-09-14 10:35:14
# @Last Modified by:   lshuns
# @Last Modified time: 2025-12-26 14:21:51

### Functions related to the correlation function 
###### calculation and plotting

import treecorr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter

def CalCorrFunc(outpath,
                cata1, cata2=None, 
                theta_Nbin=9, theta_bin_slop=0.05,
                theta_min=0.5, theta_max=300., theta_unit='arcmin', 
                theta_bin_type='Log', 
                num_threads=8, 
                c12_cata1=None, 
                c12_cata2=None,
                var_method='shot',
                Npatch=None,
                outpath_cov=None):
    """
    correlation function from treecorr
        if cata2 is None: auto correlation
        if cata2 is not None: cross correlation
    used columns:
        RA, DEC in deg
        shape_weight, e1, e2
    theta bins and slop choice:
        KV450: 9 and 0.05
        K1000: 4000 and 1.5 (Fine-binning) OR xxx and 0.08 (Broad bins)
    """

    # correct c-term and build treecorr catalogue
    if c12_cata1 is None:
        c1 = np.average(cata1['e1'].values, weights=cata1['shape_weight'].values)
        c2 = np.average(cata1['e2'].values, weights=cata1['shape_weight'].values)
        print(">>> cata1: calculated c1, c2", c1, c2)
    else:
        c1 = c12_cata1[0]
        c2 = c12_cata1[1]
        print(">>> cata1: provided c1, c2", c1, c2)
    if Npatch is None:
        cata1 = treecorr.Catalog(ra=cata1["RA"].values, dec=cata1["DEC"].values, 
                                ra_units="deg", dec_units="deg", 
                                w=cata1['shape_weight'].values,
                                g1=cata1["e1"].values - c1, 
                                g2=cata1["e2"].values - c2)
    else:
        cata1 = treecorr.Catalog(ra=cata1["RA"].values, dec=cata1["DEC"].values, 
                                ra_units="deg", dec_units="deg", 
                                w=cata1['shape_weight'].values,
                                g1=cata1["e1"].values - c1, 
                                g2=cata1["e2"].values - c2,
                                npatch=Npatch)
        print(f'+++ split sample to {Npatch} pathes with centers: {cata1.patch_centers}')
    if cata2 is not None:
        if c12_cata2 is None:
            c1 = np.average(cata2['e1'].values, weights=cata2['shape_weight'].values)
            c2 = np.average(cata2['e2'].values, weights=cata2['shape_weight'].values)
            print(">>> cata2: calculated c1, c2", c1, c2)
        else:
            c1 = c12_cata2[0]
            c2 = c12_cata2[1]
            print(">>> cata2: provided c1, c2", c1, c2)
        if Npatch is None:
            cata2 = treecorr.Catalog(ra=cata2["RA"].values, dec=cata2["DEC"].values, 
                                    ra_units="deg", dec_units="deg", 
                                    w=cata2['shape_weight'].values,
                                    g1=cata2["e1"].values - c1, 
                                    g2=cata2["e2"].values - c2)
        else:
            cata2 = treecorr.Catalog(ra=cata2["RA"].values, dec=cata2["DEC"].values, 
                                    ra_units="deg", dec_units="deg", 
                                    w=cata2['shape_weight'].values,
                                    g1=cata2["e1"].values - c1, 
                                    g2=cata2["e2"].values - c2,
                                    patch_centers=cata1.patch_centers)

    # build treecorr
    gg = treecorr.GGCorrelation(nbins=theta_Nbin, min_sep=theta_min, max_sep=theta_max, 
                                sep_units=theta_unit, 
                                bin_slop=theta_bin_slop,
                                bin_type=theta_bin_type,
                                var_method=var_method)
    # run
    gg.process(cata1, cata2, num_threads=num_threads)
    # save results
    gg.write(outpath, precision=12)
    print("TreeCorr results saved in", outpath)

    # covariance matrix
    if outpath_cov is not None:
        np.savetxt(outpath_cov, gg.cov)
        print('TreeCorr covariance saved in', outpath_cov)

def _vhlines(vORh, lines, line_styles=None, line_colors=None, line_labels=None, line_widths=None, ax=plt):
    """
    Add vertical or horizontal lines to the main plots

    Parameters
    ----------
    vORh : {'v', 'h'}
        vertical line (v) or horizontal line (h)

    lines : array-like of floats
        where to place the lines (x-axis for vline, y-axis for hline)

    line_styles : array-like of linestyles, default: 'dashed'

    line_colors : array-like of colors, default: 'k'

    line_labels : array-like of strings, default: ''

    line_widths : array-like of floats, default: 1

    ax : matplotlib Axes object, default: matplotlib.pyplot

    Returns
    -------
        None
    """

    for i, line in enumerate(lines):

        if line_styles is not None:
            line_style = line_styles[i]
        else:
            line_style = '--'

        if line_colors is not None:
            line_color = line_colors[i]
        else:
            line_color = 'k'

        if line_labels is not None:
            line_label = line_labels[i]
        else:
            line_label = ''

        if line_widths is not None:
            line_width = line_widths[i]
        else:
            line_width = 1

        if vORh == 'v':
            ax.axvline(x=line, ls=line_style, label=line_label, color=line_color, linewidth=line_width)
        elif vORh == 'h':
            ax.axhline(y=line, ls=line_style, label=line_label, color=line_color, linewidth=line_width)
        else:
            raise Exception(f'Unsupported vORh value: {vORh}')

def CorrPlotFunc_half(outpath, nzbins,
                    CATAs, which_half,
                    COLORs, 
                    yscaling = 1e4, 
                    LABELs=None, 
                    LINEs=None, LINEWs=None, 
                    POINTs=None, POINTSs=None, ERRORSIZEs=None,
                    subLABEL_locX=0.1, subLABEL_locY=0.8,
                    XRANGE=None, YRANGE=None,
                    XLABEL=None, YLABEL=None, TITLE=None,
                    xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                    vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                    hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                    xlog=False, invertX=False, ylog=False, invertY=False, 
                    loc_legend='best', legend_frame=False,
                    font_size=9, usetex=False,
                    fill_between_xs_list=None, 
                    fill_between_yLows_list=None, fill_between_yHighs_list=None,
                    fill_between_COLORs_list=None, fill_between_alphas_list=None,
                    FIGSIZE=[6.4, 4.8],
                    TIGHT=False, 
                    N_rows=None, N_cols=None,
                    timesX=True,
                    xi_p_norm=1, xi_m_norm=1,
                    xi_p_err_col='xierr_p', xi_m_err_col='xierr_m'):
    """
    plot the correlation function
        only half: which_half = 'xip' or 'xim'
    used columns in cata:
        ito, jto, theta, xi_p, xi_m, xi_p_err_col, xi_m_err_col
    """

    # font size
    plt.rc('font', size=font_size)
    # tex
    plt.rcParams["text.usetex"] = usetex

    if outpath != 'show':
        backend_orig = plt.get_backend()
        plt.switch_backend("agg")

    fig, ax = plt.subplots(nzbins, nzbins, sharex=True, sharey=True, figsize=FIGSIZE)
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    if TIGHT:
        # for x,y label
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axes
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        if TITLE is not None:
            plt.title(TITLE)

    # blank region
    for row_tmp in range(nzbins):
        for col_tmp in range(row_tmp+1, nzbins):
            ax[row_tmp, col_tmp].axis('off')

    # plot
    handles = []
    i_plot = 0
    for i in range(nzbins):
        for j in range(i, nzbins):
            row = j
            col = i

            for i_cata, cata in enumerate(CATAs):

                # get plot properties
                CR = COLORs[i_cata]
                if LABELs is not None:
                    LAB = LABELs[i_cata]
                else:
                    LAB = None
                if LINEs is not None:
                    LN = LINEs[i_cata]
                else:
                    LN = '--'
                if LINEWs is not None:
                    LW = LINEWs[i_cata]
                else:
                    LW = 1
                if POINTs is not None:
                    PI = POINTs[i_cata]
                else:
                    PI = 'o'
                if POINTSs is not None:
                    MS = POINTSs[i_cata]
                else:
                    MS = 2
                if ERRORSIZEs is not None:
                    ERRORSIZE = ERRORSIZEs[i_cata]
                else:
                    ERRORSIZE = 2

                # get values
                mask_tmp = (cata['ito']==i) & (cata['jto']==j)
                theta = cata.loc[mask_tmp, 'theta'].values
                if which_half == 'xip':
                    xi = cata.loc[mask_tmp, 'xi_p'].values / xi_p_norm
                    xierr = cata.loc[mask_tmp, xi_p_err_col].values / np.abs(xi_p_norm)
                elif which_half == 'xim':
                    xi = cata.loc[mask_tmp, 'xi_m'].values / xi_m_norm
                    xierr = cata.loc[mask_tmp, xi_m_err_col].values / np.abs(xi_m_norm)
                else:
                    raise Exception(f'Unsupported which_half {which_half}')
                del mask_tmp, cata

                # plot
                if timesX:
                    tmp = ax[row, col].errorbar(theta, theta*xi*yscaling, yerr=np.abs(theta*yscaling)*xierr,
                                                        color=CR, 
                                                        label=LAB, 
                                                        linestyle=LN, linewidth=LW, 
                                                        marker=PI, markersize=MS, capsize=ERRORSIZE)
                else:
                    tmp = ax[row, col].errorbar(theta, xi*yscaling, yerr=xierr*np.abs(yscaling),
                                                        color=CR, 
                                                        label=LAB, 
                                                        linestyle=LN, linewidth=LW, 
                                                        marker=PI, markersize=MS, capsize=ERRORSIZE)
                if (i_plot==0):
                    handles.append(tmp[0])
                del tmp

            if fill_between_xs_list is not None:
                fill_between_xs = fill_between_xs_list[i_plot]
                fill_between_yLows = fill_between_yLows_list[i_plot]
                fill_between_yHighs = fill_between_yHighs_list[i_plot]
                fill_between_COLORs = fill_between_COLORs_list[i_plot]
                fill_between_alphas = fill_between_alphas_list[i_plot]
                for i_fill, fill_between_x in enumerate(fill_between_xs):
                    ax[row, col].fill_between(fill_between_x, 
                                fill_between_yLows[i_fill], fill_between_yHighs[i_fill],
                                alpha=fill_between_alphas[i_fill],
                                color=fill_between_COLORs[i_fill])
            i_plot += 1

            LABEL = f'{i+1}-{j+1}'
            ax[row, col].text(subLABEL_locX, subLABEL_locY, LABEL, transform=ax[row, col].transAxes)

            if XRANGE is not None:
                ax[row, col].set_xlim(XRANGE[0], XRANGE[1])
            if YRANGE is not None:
                ax[row, col].set_ylim(YRANGE[0], YRANGE[1])

            if vlines is not None:
                _vhlines('v', vlines, line_styles=vline_styles, line_colors=vline_colors, line_labels=vline_labels, line_widths=vline_widths, ax=ax[row, col])
            if hlines is not None:
                _vhlines('h', hlines, line_styles=hline_styles, line_colors=hline_colors, line_labels=hline_labels, line_widths=hline_widths, ax=ax[row, col])

            if xlog:
                ax[row, col].set_xscale('log')
                ax[row, col].xaxis.set_minor_formatter(NullFormatter())
            if ylog:
                ax[row, col].set_yscale('log')
                ax[row, col].yaxis.set_minor_formatter(NullFormatter())

            if xtick_min_label:
                if xlog:
                    ax[row, col].xaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))
                else:
                    ax[row, col].xaxis.set_minor_locator(AutoMinorLocator())
            if ytick_min_label:
                if ylog:
                    ax[row, col].yaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))
                else:
                    ax[row, col].yaxis.set_minor_locator(AutoMinorLocator())

            if xtick_spe is not None:
                ax[row, col].set_xticks(xtick_spe[0])
                ax[row, col].set_xticklabels(xtick_spe[1])
            if ytick_spe is not None:
                ax[row, col].set_yticks(ytick_spe[0])
                ax[row, col].set_yticklabels(ytick_spe[1])

            if invertY:
                plt.gca().invert_yaxis()
            if invertX:
                plt.gca().invert_xaxis()

    if not TIGHT:
        fig.text(0.5, 0.03, XLABEL, ha='center')

        fig.text(0.07, 0.5, YLABEL, 
             horizontalalignment='center', verticalalignment='center',
             rotation='vertical')

        if TITLE is not None:
            fig.text(0.5, 0.90, TITLE, ha='center')

    if LABELs is not None:
        fig.legend(handles, LABELs, 
                loc = 'upper right',
                bbox_to_anchor=(0.92, 0.92), fancybox=False, shadow=False, frameon=False)

    if TIGHT:
        plt.tight_layout()

    if outpath == 'show':
        plt.show()
        plt.close()
    else:
        plt.savefig(outpath, dpi=300)
        plt.close()
        plt.switch_backend(backend_orig)
        print("Errorbar plot saved as", outpath)

def CorrPlotFunc_half_diff2base(outpath, nzbins,
                    CATAs, CATA_base, which_half,
                    COLORs, 
                    yscaling = 1e4, 
                    LABELs=None, 
                    LINEs=None, LINEWs=None, 
                    POINTs=None, POINTSs=None, ERRORSIZEs=None,
                    subLABEL_locX=0.1, subLABEL_locY=0.8,
                    XRANGE=None, YRANGE=None,
                    XLABEL=None, YLABEL=None, TITLE=None,
                    xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                    vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                    hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                    xlog=False, invertX=False, ylog=False, invertY=False, 
                    loc_legend='best', legend_frame=False,
                    font_size=9, usetex=False,
                    fill_between_xs_list=None, 
                    fill_between_yLows_list=None, fill_between_yHighs_list=None,
                    fill_between_COLORs_list=None, fill_between_alphas_list=None,
                    FIGSIZE=[6.4, 4.8],
                    TIGHT=False, 
                    N_rows=None, N_cols=None,
                    timesX=True,
                    xi_p_norm=1, xi_m_norm=1,
                    xi_p_err_col='xierr_p', xi_m_err_col='xierr_m'):
    """
    plot the correlation function difference with respect to a given fiducial
        only half: which_half = 'xip' or 'xim'
    The difference between 
    used columns in cata:
        ito, jto, theta, xi_p, xi_m, xi_p_err_col, xi_m_err_col

    """

    # font size
    plt.rc('font', size=font_size)
    # tex
    plt.rcParams["text.usetex"] = usetex

    if outpath != 'show':
        backend_orig = plt.get_backend()
        plt.switch_backend("agg")

    fig, ax = plt.subplots(nzbins, nzbins, sharex=True, sharey=True, figsize=FIGSIZE)
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    if TIGHT:
        # for x,y label
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axes
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.grid(False)
        plt.xlabel(XLABEL)
        plt.ylabel(YLABEL)
        if TITLE is not None:
            plt.title(TITLE)

    # blank region
    for row_tmp in range(nzbins):
        for col_tmp in range(row_tmp+1, nzbins):
            ax[row_tmp, col_tmp].axis('off')

    # plot
    handles = []
    i_plot = 0
    for i in range(nzbins):
        for j in range(i, nzbins):
            row = j
            col = i

            for i_cata, cata in enumerate(CATAs):

                # get plot properties
                CR = COLORs[i_cata]
                if LABELs is not None:
                    LAB = LABELs[i_cata]
                else:
                    LAB = None
                if LINEs is not None:
                    LN = LINEs[i_cata]
                else:
                    LN = '--'
                if LINEWs is not None:
                    LW = LINEWs[i_cata]
                else:
                    LW = 1
                if POINTs is not None:
                    PI = POINTs[i_cata]
                else:
                    PI = 'o'
                if POINTSs is not None:
                    MS = POINTSs[i_cata]
                else:
                    MS = 2
                if ERRORSIZEs is not None:
                    ERRORSIZE = ERRORSIZEs[i_cata]
                else:
                    ERRORSIZE = 2

                # get values from cata
                mask_tmp = (cata['ito']==i) & (cata['jto']==j)
                theta = cata.loc[mask_tmp, 'theta'].values
                if which_half == 'xip':
                    xi = cata.loc[mask_tmp, 'xi_p'].values / xi_p_norm
                    xierr = cata.loc[mask_tmp, xi_p_err_col].values / np.abs(xi_p_norm)
                elif which_half == 'xim':
                    xi = cata.loc[mask_tmp, 'xi_m'].values / xi_m_norm
                    xierr = cata.loc[mask_tmp, xi_m_err_col].values / np.abs(xi_m_norm)
                else:
                    raise Exception(f'Unsupported which_half {which_half}')
                del mask_tmp, cata

                # get values from base
                mask_tmp = (CATA_base['ito']==i) & (CATA_base['jto']==j)
                theta_base = CATA_base.loc[mask_tmp, 'theta'].values
                if which_half == 'xip':
                    xi_base = CATA_base.loc[mask_tmp, 'xi_p'].values / xi_p_norm
                elif which_half == 'xim':
                    xi_base = CATA_base.loc[mask_tmp, 'xi_m'].values / xi_m_norm
                else:
                    raise Exception(f'Unsupported which_half {which_half}')
                del mask_tmp

                # calculate the difference
                theta_plot = (theta_base + theta)/2.
                xi_plot = xi - xi_base
                del theta_base, theta, xi_base, xi

                # plot
                if timesX:
                    tmp = ax[row, col].errorbar(theta_plot, theta_plot*xi_plot*yscaling, 
                                                yerr=np.abs(theta_plot*yscaling)*xierr,
                                                        color=CR, 
                                                        label=LAB, 
                                                        linestyle=LN, linewidth=LW, 
                                                        marker=PI, markersize=MS, capsize=ERRORSIZE)
                else:
                    tmp = ax[row, col].errorbar(theta_plot, xi_plot*yscaling, 
                                                yerr=xierr*np.abs(yscaling),
                                                        color=CR, 
                                                        label=LAB, 
                                                        linestyle=LN, linewidth=LW, 
                                                        marker=PI, markersize=MS, capsize=ERRORSIZE)
                if (i_plot==0):
                    handles.append(tmp[0])
                del tmp

            if fill_between_xs_list is not None:
                fill_between_xs = fill_between_xs_list[i_plot]
                fill_between_yLows = fill_between_yLows_list[i_plot]
                fill_between_yHighs = fill_between_yHighs_list[i_plot]
                fill_between_COLORs = fill_between_COLORs_list[i_plot]
                fill_between_alphas = fill_between_alphas_list[i_plot]
                for i_fill, fill_between_x in enumerate(fill_between_xs):
                    ax[row, col].fill_between(fill_between_x, 
                                fill_between_yLows[i_fill], fill_between_yHighs[i_fill],
                                alpha=fill_between_alphas[i_fill],
                                color=fill_between_COLORs[i_fill])
            i_plot += 1

            LABEL = f'{i+1}-{j+1}'
            ax[row, col].text(subLABEL_locX, subLABEL_locY, LABEL, transform=ax[row, col].transAxes)

            if XRANGE is not None:
                ax[row, col].set_xlim(XRANGE[0], XRANGE[1])
            if YRANGE is not None:
                ax[row, col].set_ylim(YRANGE[0], YRANGE[1])

            if vlines is not None:
                _vhlines('v', vlines, line_styles=vline_styles, line_colors=vline_colors, line_labels=vline_labels, line_widths=vline_widths, ax=ax[row, col])
            if hlines is not None:
                _vhlines('h', hlines, line_styles=hline_styles, line_colors=hline_colors, line_labels=hline_labels, line_widths=hline_widths, ax=ax[row, col])

            if xlog:
                ax[row, col].set_xscale('log')
                ax[row, col].xaxis.set_minor_formatter(NullFormatter())
            if ylog:
                ax[row, col].set_yscale('log')
                ax[row, col].yaxis.set_minor_formatter(NullFormatter())

            if xtick_min_label:
                if xlog:
                    ax[row, col].xaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))
                else:
                    ax[row, col].xaxis.set_minor_locator(AutoMinorLocator())
            if ytick_min_label:
                if ylog:
                    ax[row, col].yaxis.set_minor_locator(LogLocator(base=10.0, subs=None, numticks=10))
                else:
                    ax[row, col].yaxis.set_minor_locator(AutoMinorLocator())

            if xtick_spe is not None:
                ax[row, col].set_xticks(xtick_spe[0])
                ax[row, col].set_xticklabels(xtick_spe[1])
            if ytick_spe is not None:
                ax[row, col].set_yticks(ytick_spe[0])
                ax[row, col].set_yticklabels(ytick_spe[1])

            if invertY:
                plt.gca().invert_yaxis()
            if invertX:
                plt.gca().invert_xaxis()

    if not TIGHT:
        fig.text(0.5, 0.03, XLABEL, ha='center')

        fig.text(0.07, 0.5, YLABEL, 
             horizontalalignment='center', verticalalignment='center',
             rotation='vertical')

        if TITLE is not None:
            fig.text(0.5, 0.90, TITLE, ha='center')

    if LABELs is not None:
        fig.legend(handles, LABELs, 
                loc = 'upper right',
                bbox_to_anchor=(0.92, 0.92), fancybox=False, shadow=False, frameon=False)

    if TIGHT:
        plt.tight_layout()

    if outpath == 'show':
        plt.show()
        plt.close()
    else:
        plt.savefig(outpath, dpi=300)
        plt.close()
        plt.switch_backend(backend_orig)
        print("Errorbar plot saved as", outpath)