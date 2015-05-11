#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
# sns.set_context('paper')
sns.set_context('poster')
from glob import glob


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
    dfB = pd.DataFrame(b, columns=['w', 'Excitation potential', 'log gf', 'EW', 'Abundance'])


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
        sns.regplot(x[x==4.5], y[x==4.5], x_estimator=np.median, fit_reg=False)
        # EP cut
        sns.regplot(x[x!=4.5], y[x!=4.5], x_estimator=np.median, truncate=True,
                label='EP cut')
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

    # plt.savefig('figures/solar_parameters_10runs.pdf')
    plt.show()


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


    c = sns.color_palette()
    f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3, xticks=np.arange(5, 6.2, 0.1),
             xticklabels=list(np.arange(5, 6.1, 0.1)) + ['No cut'])

    sns.regplot('EP', 'Teff', df, truncate=True, order=2, ax=ax1).set_xlabel('')
    ax1.plot(6.1, data[-1, 1], 'o', color=c[0], ms=6)
    ax1.set_ylabel(r'$\mathrm{T_{eff}}$ [K]')

    sns.regplot('EP', 'feh', df, truncate=True, order=2, ax=ax2).set_xlabel('')
    ax2.plot(6.1, data[-1, 4], 'o', color=c[1], ms=6)
    ax2.set_ylabel('[Fe/H]')

    sns.regplot('EP', 'vt', df, truncate=True, order=2, ax=ax3)
    ax3.plot(6.1, data[-1, 3], 'o', color=c[2], ms=6)
    ax3.set_ylabel(r'$\xi_\mathrm{micro}$ [km/s]')

    plt.axis((4.95, 6.15, -0.5, 3.3))
    plt.show()
    # plt.savefig('figures/HD20010_parameters_cuts.pdf')


def main():
    """Main function
    :returns: TODO
    """
    # fig_abundance()
    # fig_EPcut_sun()
    fig_HD20010_parameters()


if __name__ == '__main__':
    main()
