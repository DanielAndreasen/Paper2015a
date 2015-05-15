#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
sns.set_style('white')
# sns.set_context('paper', font_scale=1.5)
sns.set_context('poster')
from glob import glob
from astropy.io import fits
from astropy import constants as c
from astropy import units as u


def fig_abundance(fout=None):
    """Figures of abundance before and after recalibration

    :fname: Name of input file
    :mean: mean value of abundance
    :fout: file name the figures should be written to
    """

    p = '../../../programs/moog/results/'
    fouts = ('figures/EWvsEP', 'figures/EWvsEP_cut')
    b = np.loadtxt(p + 'Fe1_PreSynth_rec.log', skiprows=1)
    b = np.c_[b[:, 0], b[:, 1], b[:, 2], b[:, 3], b[:, 5]]
    dfB = pd.DataFrame(b, columns=['w', 'Excitation potential',
                                   'log gf', 'EW', 'Abundance'])

    idx = dfB.Abundance < 8.47
    sns.interactplot('Excitation potential', 'EW', 'Abundance', dfB,
                     filled=True, levels=100,
                     scatter_kws={'ms': 0})

    sns.interactplot('Excitation potential', 'EW', 'Abundance', dfB[idx],
                     filled=True, levels=100, colorbar=False,
                     scatter_kws={'label': 'Good abundance'})

    sns.interactplot('Excitation potential', 'EW', 'Abundance', dfB[~idx],
                     filled=True, levels=100, colorbar=False,
                     scatter_kws={'alpha': 0.3, 'ms': 4,
                                  'label': 'Bad abundance',
                                  'color': 'r'})
    plt.legend(loc='best')
    plt.savefig('%s.pdf' % fouts[0], format='pdf')

    a = np.loadtxt(p + 'Fe1_PostSynth_cut.log', skiprows=1)
    ax = sns.jointplot(a[:, 1], a[:, 3], stat_func=None, kind='kde')
    ax.set_axis_labels(xlabel='Excitation potential',
                       ylabel='\nEW')
    # plt.show()
    plt.savefig('%s.pdf' % fouts[1], format='pdf')


def fig_EPcut_sun(fout=None):
    """Figure for the parameters of the Sun as function of the cut in EP

    :fname: file name with result
    :fout: File name figure is written to
    """

    def _get_params(fname):
        with open(fname, 'r') as lines:
            for _ in range(2):
                lines.readline()
            line = lines.readline().strip().split()
            if line[0] == '#':
                line = lines.readline().strip().split()
        T = int(line[1])
        logg = float(line[4])
        vt = float(line[6])
        if len(line) == 8:
            feh = float(line[-1][4:])
        else:
            feh = float(line[-1])
        return T, logg, vt, feh

    def _plot_result(data, xlabel=None, ylabel=None, solar=None):
        x, y = data
        # No EP cut
        sns.regplot(x[x == 4.5], y[x == 4.5], x_estimator=np.median,
                    fit_reg=False)
        # EP cut
        sns.regplot(x[x != 4.5], y[x != 4.5], x_estimator=np.median,
                    truncate=True, label='EP cut')
        plt.xticks([4.5, 5.0, 5.5], ['No cut', '5.0', '5.5'])
        if ylabel:
            plt.ylabel(ylabel)
        if solar:
            plt.hlines(solar, 4.5, 5.5)

    p = '/home/daniel/Software/SPECPAR/Sun/'
    files = glob(p + 'Out_moog_*')
    parameters = np.zeros((len(files) - 1, 5))

    i = 0
    for file in files:
        if 'Poisson' not in file:
            solar = _get_params(file)
        else:
            if 'filtered' in file:
                epcut = float(file.split('filtered_2_')[1].strip('.moog'))
            else:
                epcut = 4.5
            T, logg, vt, feh = _get_params(file)
            parameters[i, 0] = T
            parameters[i, 1] = logg
            parameters[i, 2] = vt
            parameters[i, 3] = feh
            parameters[i, 4] = epcut
            i += 1

    T, logg, vt, feh, epcut = parameters.T
    ax1 = plt.subplot(221)
    _plot_result(data=(epcut, T), ylabel=r'$\mathrm{T_{eff}}$', solar=solar[0])

    ax2 = plt.subplot(222, sharex=ax1)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    _plot_result(data=(epcut, logg), ylabel=r'$\log(g)$', solar=solar[1])

    ax3 = plt.subplot(223, sharex=ax1)
    _plot_result(data=(epcut, vt), xlabel=True, ylabel=r'$\xi_\mathrm{micro}$',
                 solar=solar[2])

    ax4 = plt.subplot(224, sharex=ax1)
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position('right')
    _plot_result(data=(epcut, feh), xlabel=True, ylabel='[Fe/H]')
    plt.hlines(0.0, 4.5, 5.5)

    plt.savefig('figures/solar_parameters_10runs.pdf')
    # plt.show()


def fig_HD20010_parameters():
    """
    Derived parameters for HD20010 at different cuts in EP
    """
    def _get_parameters(file):
        with open(file, 'r') as lines:
            for _ in range(2):
                lines.readline()
            line = lines.readline()
        z = line.strip().replace(' ', '').split('=')
        return z

    p = '/home/daniel/Software/SPECPAR/HD20010/'
    N = len(glob(p+'Out*.moog'))
    data = np.zeros((N, 5))

    prefixes = (p+'Out_moog_b.', p+'Out_moog_')
    with open(p+'star_file', 'r') as files:
        # Dealing with the header.
        for _ in range(2):
            files.readline()
        for index, file in enumerate(files):
            ep = file.split('.moog')[0].split('_')[-1]
            ep = float(ep)
            file = file.split(' ')[0]
            file1 = prefixes[0] + file
            file2 = prefixes[1] + file
            if os.path.isfile(file1):
                z = _get_parameters(file1)
            else:
                z = _get_parameters(file2)
            try:
                if ep == 1.0:
                    ep = 6.2
                data[index, 0] = ep
                data[index, 1] = int(z[1][:4])  # Teff
                data[index, 2] = float(z[2][:4])  # logg
                data[index, 3] = float(z[3][:4])  # microturbulence
                data[index, 4] = float(z[-1])  # metalicity
            except ValueError:
                data[index, 0] = ep
                data[index, 1] = int(z[1][:4])  # Teff
                data[index, 2] = float(z[2][:4])  # logg
                data[index, 3] = np.nan  # microturbulence
                data[index, 4] = float(z[-1])  # metalicity

    df = pd.DataFrame({'EP':   data[0:-1, 0],
                       'Teff': data[0:-1, 1],
                       'logg': data[0:-1, 2],
                       'vt':   data[0:-1, 3],
                       'feh':  data[0:-1, 4]})

    xlim1 = [4.95, 6.05]
    xlim2 = [6.95, 7.05]
    xlim1ratio = (xlim1[1]-xlim1[0])/(xlim2[1]-xlim2[0]+xlim1[1]-xlim1[0])
    xlim2ratio = (xlim2[1]-xlim2[0])/(xlim2[1]-xlim2[0]+xlim1[1]-xlim1[0])
    gs = gridspec.GridSpec(3, 2, width_ratios=[xlim1ratio, xlim2ratio])

    c = sns.color_palette()

    fig = plt.figure()
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax1.plot(xlim1, [6160]*2, '--k')
    sns.regplot('EP', 'Teff', df, truncate=True, order=2, ax=ax1,
                color=c[0]).set_xlabel('')
    ax1.set_ylabel('Teff [K]')
    ax2.plot(7.0, data[-1, 1], 'o', color=c[0], ms=6)
    ax1.set_xlim(xlim1)
    ax2.set_xlim(xlim2)
    plt.subplots_adjust(wspace=0.03)

    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright='off')
    plt.setp(ax2, xticks=[7.0], xticklabels=['No EP cut'])

    ax2.xaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)
    kwargs = dict(color='k', clip_on=False, linewidth=1)

    ylim = (5800, 7400)
    dx = 0.01 * (xlim1[1]-xlim1[0])/xlim1ratio
    dy = 0.02 * (ylim[1]-ylim[0])
    ax1.plot((xlim1[1]-dx, xlim1[1]+dx), (ylim[1]-dy, ylim[1]+dy), **kwargs)
    ax1.plot((xlim1[1]-dx, xlim1[1]+dx), (ylim[0]-dy, ylim[0]+dy), **kwargs)
    ax2.plot((xlim2[0]-dx, xlim2[0]+dx), (ylim[1]-dy, ylim[1]+dy), **kwargs)
    ax2.plot((xlim2[0]-dx, xlim2[0]+dx), (ylim[0]-dy, ylim[0]+dy), **kwargs)

    ax1.set_ylim(ylim)
    ax2.set_ylim(ylim)

    # [Fe/H]
    ax3 = fig.add_subplot(gs[2])
    ax4 = fig.add_subplot(gs[3])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)

    ax3.plot(xlim1, [-0.23]*2, '--k')
    sns.regplot('EP', 'feh', df, truncate=True, order=2, ax=ax3,
                color=c[1]).set_xlabel('')
    ax3.set_ylabel('[Fe/H]')
    ax4.plot(7.0, data[-1, 4], 'o', color=c[1], ms=6)
    ax3.set_xlim(xlim1)
    ax4.set_xlim(xlim2)
    plt.subplots_adjust(wspace=0.03)

    ax3.spines['right'].set_visible(False)
    ax4.spines['left'].set_visible(False)
    ax4.yaxis.tick_right()
    ax4.tick_params(labelright='off')
    plt.setp(ax4, xticks=[7.0], xticklabels=['No EP cut'])

    ax4.xaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)
    kwargs = dict(color='k', clip_on=False, linewidth=1)

    ylim = (-0.30, 0.00)
    dx = 0.01 * (xlim1[1]-xlim1[0])/xlim1ratio
    dy = 0.02 * (ylim[1]-ylim[0])
    ax3.plot((xlim1[1]-dx, xlim1[1]+dx), (ylim[1]-dy, ylim[1]+dy), **kwargs)
    ax3.plot((xlim1[1]-dx, xlim1[1]+dx), (ylim[0]-dy, ylim[0]+dy), **kwargs)
    ax4.plot((xlim2[0]-dx, xlim2[0]+dx), (ylim[1]-dy, ylim[1]+dy), **kwargs)
    ax4.plot((xlim2[0]-dx, xlim2[0]+dx), (ylim[0]-dy, ylim[0]+dy), **kwargs)

    ax3.set_ylim(ylim)
    ax4.set_ylim(ylim)

    # micorturbulence
    ax5 = fig.add_subplot(gs[4])
    ax6 = fig.add_subplot(gs[5])
    # plt.setp(ax5.get_xticklabels(), visible=False)
    # plt.setp(ax6.get_xticklabels(), visible=False)

    sns.regplot('EP', 'vt', df, truncate=True, order=2, ax=ax5,
                color=c[2]).set_xlabel('')
    ax5.set_ylabel(r'$\xi_\mathrm{micro}$ [km/s]')
    ax6.plot(7.0, data[-1, 3], 'o', color=c[2], ms=6)
    ax5.set_xlim(xlim1)
    ax6.set_xlim(xlim2)
    plt.subplots_adjust(wspace=0.03)

    ax5.spines['right'].set_visible(False)
    ax6.spines['left'].set_visible(False)
    ax6.yaxis.tick_right()
    ax6.tick_params(labelright='off')
    plt.setp(ax6, xticks=[7.0], xticklabels=['No EP cut'])

    ax6.xaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)
    kwargs = dict(color='k', clip_on=False, linewidth=1)

    ylim = (0.20, 3.10)
    dx = 0.01 * (xlim1[1]-xlim1[0])/xlim1ratio
    dy = 0.02 * (ylim[1]-ylim[0])
    ax5.plot((xlim1[1]-dx, xlim1[1]+dx), (ylim[1]-dy, ylim[1]+dy), **kwargs)
    ax5.plot((xlim1[1]-dx, xlim1[1]+dx), (ylim[0]-dy, ylim[0]+dy), **kwargs)
    ax6.plot((xlim2[0]-dx, xlim2[0]+dx), (ylim[1]-dy, ylim[1]+dy), **kwargs)
    ax6.plot((xlim2[0]-dx, xlim2[0]+dx), (ylim[0]-dy, ylim[0]+dy), **kwargs)

    ax5.set_ylim(ylim)
    ax6.set_ylim(ylim)

    # plt.show()
    plt.savefig('figures/HD20010_parameters_cuts.pdf')


def fig_spectral_region():
    """
    Plot a high Teff synthetic spectra and a low Teff synthetic spectra
    """

    path = '/home/daniel/Documents/Uni/phdproject/data/HD20010/article/figures/'
    suffix = '.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
    w = fits.getdata('%sWAVE_PHOENIX-ACES-AGSS-COND-2011.fits' % path)
    fl = fits.getdata('%slte02700-4.50-0.0%s' % (path, suffix))
    fm = fits.getdata('%slte03500-4.00-0.0%s' % (path, suffix))
    fh = fits.getdata('%slte06200-4.00-0.0%s' % (path, suffix))
    h = fits.getheader('%slte02700-4.50-0.0%s' % (path, suffix))

    r = 500
    w = w[::r]
    fl = fl[::r]
    fm = fm[::r]
    fh = fh[::r]

    from astropy.analytic_functions import blackbody_lambda
    _wmax = lambda Teff: (c.b_wien/Teff).to('AA')
    PHOENIX_unit = u.erg/(u.s * (u.cm**2) * u.cm * u.sr)
    # PHOENIX_unit = u.erg/(u.s * (u.cm**2) * u.cm)

    Teff = np.arange(2000, 8500, 50) * u.K
    wmax = _wmax(Teff)
    bb_max = [blackbody_lambda(wi, Ti) for wi, Ti in zip(wmax, Teff)]
    Bm = np.array([bi.value * 4*np.pi for bi in bb_max]) * u.erg/(u.s * (u.cm**2) * u.AA)
    # Bm = np.array([bi.value for bi in bb_max]) * bb_max[0].unit
    # Bm_new = Bm.to(PHOENIX_unit)
    # Bm_new = Bm.to(PHOENIX_unit, equivalencies=u.spectral_density(wmax))

    print 'PHOENIX unit: %s' % h['BUNIT']
    print 'Astropy unit: %s' % Bm.unit
    # print '    New unit: %s' % Bm_new.unit


    plt.plot(w, fl, label=r'$T_\mathrm{eff}=2700\mathrm{K, }\log g=4.50$')
    # plt.plot(w, blackbody_lambda(w, 2700*u.K).to(PHOENIX_unit))
    plt.plot(w, fm, label=r'$T_\mathrm{eff}=3500\mathrm{K, }\log g=4.00$')
    # plt.plot(w, blackbody_lambda(w, 3500*u.K).to(PHOENIX_unit))
    plt.plot(w, fh, label=r'$T_\mathrm{eff}=6200\mathrm{K, }\log g=4.00$')
    # plt.plot(w, blackbody_lambda(w, 6200*u.K).to(PHOENIX_unit) * 4.5)
    plt.plot(wmax.value, Bm.value, '-r')
    # plt.plot(wmax.value, Bm_new.value, '-k')

    plt.xlabel(r'$\lambda$ ({0})'.format(wmax.unit))
    # plt.ylabel(r'Flux [{0}]'.format(Bm_new.unit))
    plt.legend(loc='best', frameon=False)
    plt.show()







def main():
    """Main function
    :returns: TODO
    """
    # fig_abundance()
    # fig_EPcut_sun()
    # fig_HD20010_parameters()
    fig_spectral_region()


if __name__ == '__main__':
    main()
