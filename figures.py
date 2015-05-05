#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_context('paper')
# sns.set_context('poster')
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
    ax = sns.jointplot(a[:, 1], a[:, 3], stat_func=None)
    ax.set_axis_labels(xlabel='Excitation potential',
                       ylabel='\nEW')
    plt.savefig('%s.pdf' % fouts[1], format='pdf')


    # plt.show()


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

    for i, file in enumerate(files):
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

    T, logg, vt, feh, epcut = parameters.T
    ax1 = plt.subplot(221)
    _plot_result(data=(epcut, T), ylabel='Teff', solar=solar[0])

    ax2 = plt.subplot(222, sharex=ax1)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    _plot_result(data=(epcut, logg), ylabel='logg', solar=solar[1])

    ax3 = plt.subplot(223, sharex=ax1)
    _plot_result(data=(epcut, vt), xlabel=True, ylabel='vt', solar=solar[2])

    ax4 = plt.subplot(224, sharex=ax1)
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position('right')
    _plot_result(data=(epcut, feh), xlabel=True, ylabel='[Fe/H]')
    plt.hlines(0.0, 4.5, 5.5)

    # plt.savefig('figures/solar_parameters_10runs.pdf')
    plt.show()


def fig_something(fout=None):
    """Figure of something

    :fout: Save file to this
    """


def main():
    """Main function
    :returns: TODO
    """
    # fig_abundance()
    fig_EPcut_sun()


if __name__ == '__main__':
    main()
