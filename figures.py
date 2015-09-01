#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib
import seaborn as sns
sns.set_style('white')
# sns.set_context('paper', font_scale=1.5)
# sns.set_context('poster')
sns.set_context('talk', font_scale=1.2)
from glob import glob
from astropy.io import fits
from astropy import constants as c
from astropy import units as u
from plot_fits import get_wavelength


matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)
params = {'legend.fontsize': 24}
matplotlib.rcParams.update(params)
matplotlib.rcParams.update({'font.size': 17})


def fig_EWvsEP(fout=None):
    """Figures of EW vs EP"""

    p = '../../../programs/moog/results/'

    df1 = pd.read_csv('figures/Fe1.dat', delimiter=r'\s+')
    df1.rename(columns={'abund': 'Abundance'}, inplace=True)
    df2 = pd.read_csv('figures/Fe2.dat', delimiter=r'\s+')
    df2.rename(columns={'abund': 'Abundance'}, inplace=True)

    g = sns.jointplot('EP', 'EW', df1, stat_func=None, kind='scatter', space=0)
    # stackoverflow.com/questions/31539815/plotting-two-distributions-in-seaborn-jointplot
    g.x = df2.EP
    g.y = df2.EW
    g.plot_joint(plt.scatter, c='r', s=30, alpha=0.6)
    g.set_axis_labels(xlabel='Excitation potential [eV]', ylabel=r'EW [m$\AA$]')
    plt.ylim(-5, 210)
    plt.show()

    plt.savefig('figures/EWvsEP.pdf', format='pdf')


def fig_solarspectrum():
    """Figure of the solar spectrum with the distribution of iron lines on top"""
    color = sns.color_palette()
    p_spec = '/home/daniel/Documents/Uni/phdproject/data/atlas/BASS2000/solarspectrum_01.fits'
    p = '../../../programs/moog/results/'
    df1 = pd.read_csv('figures/Fe1.dat', delimiter=r'\s+')
    df1.rename(columns={'abund': 'Abundance'}, inplace=True)
    df2 = pd.read_csv('figures/Fe2.dat', delimiter=r'\s+')
    df2.rename(columns={'abund': 'Abundance'}, inplace=True)

    I = fits.getdata(p_spec)
    I /= np.median(I)
    I *= 40
    w = get_wavelength(fits.getheader(p_spec))
    w, I = w[(w>9500) & (w<25000)], I[(w>9500) & (w<25000)]
    plt.plot(w[::50], I[::50], '-k', alpha=0.3)
    plt.hist(df1.wavelength, color=color[1], bins=30, label='FeI')
    plt.hist(df2.wavelength, color=color[0], bins=20, label='FeII')
    plt.xlabel(r'Wavelength $\AA$', fontsize=24)
    plt.ylabel('Number of lines', fontsize=24)
    plt.title('Recalibrated iron lines', fontsize=27)
    plt.legend(loc=2, frameon=False)
    fig = plt.gcf()
    # fig.subplots_adjust(bottom=0.15)
    plt.xticks(np.arange(8000, 26001, 3000), np.arange(8000, 26001, 3000))
    plt.tight_layout()
    plt.savefig('figures/EWvsEP_cut.pdf', format='pdf')
    # plt.show()


def fig_solarparams():
    '''Solar parameters at different SNR'''
    df = pd.read_csv('solar_snr_params.dat', delimiter=r'\s+')
    snrs = sorted(list(set(df.SNR)))
    params = np.zeros((len(snrs), 9))

    fig = plt.figure()
    for i, snr in enumerate(snrs):
        # Teff
        t = np.mean(df['Teff'][df.SNR == snr])
        s = 3*np.std(df['Teff'][df.SNR == snr])
        ax1 = fig.add_subplot(221)
        ax1.errorbar(snr, t-5777, yerr=s, fmt='ok')
        ax1.set_xticklabels([])
        params[i, 0] = snr
        params[i, 1] = t
        params[i, 2] = s

        # logg
        t = np.mean(df['logg'][df.SNR == snr])
        s = 3*np.std(df['logg'][df.SNR == snr])
        ax2 = fig.add_subplot(222)
        plt.errorbar(snr, t-4.438, yerr=s, fmt='ok')
        ax2.set_xticklabels([])
        params[i, 3] = t
        params[i, 4] = s

        # feh
        t = np.mean(df['feh'][df.SNR == snr])
        s = 3*np.std(df['feh'][df.SNR == snr])
        ax3 = fig.add_subplot(223)
        ax3.errorbar(snr, t, yerr=s, fmt='ok')
        params[i, 5] = t
        params[i, 6] = s

        # vt
        t = np.mean(df['vt'][df.SNR == snr])
        s = 3*np.std(df['vt'][df.SNR == snr])
        ax4 = fig.add_subplot(224)
        ax4.errorbar(snr, t-1.0, yerr=s, fmt='ok')
        params[i, 7] = t
        params[i, 8] = s

    # for i in range(params.shape[0]):
    #     print "{0:.0f} &  ${1:.0f} \\pm {2:.0f}$  & ${3:.2f} \\pm {4:.2f}$ & ${5:.2f} \\pm {6:.2f}$ & ${7:.2f} \\pm {8:.2f}$ \\\\".format(*params[i, :])

    ax1.grid(True)
    ax1.set_ylabel('Teff [K] - 5777K')
    ax1.set_xlim((15, 320))
    ax1.set_ylim(-100, 200)
    ax1.set_xticks([25, 50, 100, 150, 225, 300])

    ax2.grid(True)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.set_ylabel('logg - 4.438')
    ax2.set_xlim((15, 320))
    ax2.set_ylim(-0.3, 0.3)
    ax2.set_xticks([25, 50, 100, 150, 225, 300])

    ax3.grid(True)
    ax3.set_ylabel('[Fe/H]')
    ax3.set_xlabel('SNR')
    ax3.set_xlim((15, 320))
    ax3.set_ylim(-0.1, 0.1)
    ax3.set_xticks([25, 50, 100, 150, 225, 300])

    ax4.grid(True)
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position('right')
    ax4.set_ylabel(r'$\xi_\mathrm{micro}$ [km/s] - 1.0km/s')
    ax4.set_xlabel('SNR')
    ax4.set_xlim((15, 320))
    ax4.set_ylim(-1, 1)
    ax4.set_xticks([25, 50, 100, 150, 225, 300])


    plt.tight_layout()
    plt.savefig('figures/solar_parameters_snr.pdf')
    # plt.show()


def fig_HD20010_parameters():
    """
    Derived parameters for HD20010 at different cuts in EW
    """

    df = pd.read_csv('HD20010_params_EWcut_fixlogg_3sigma.dat')
    c = sns.color_palette()

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.errorbar(df.EW, df.Teff, yerr=df.Tefferr, fmt='o', color=c[0])
    ax1.fill_between([0, 20], [6131-255]*2, [6131+255]*2, color='k', alpha=0.2, edgecolor='w')
    ax1.hlines(6131, 0, 20, linestyle='--')
    ax1.set_xlim(-1, 21)
    ax1.set_ylim(5800, 7800)
    ax1.set_ylabel('Teff [K]')

    ax2 = fig.add_subplot(312)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.errorbar(df.EW, df.feh, yerr=df.feherr, fmt='o', color=c[1])
    ax2.fill_between([0, 20], [-0.23-0.14]*2, [-0.23+0.14]*2, color='k', alpha=0.2, edgecolor='w')
    ax2.hlines(-0.23, 0, 20, linestyles='--')
    ax2.set_xlim(-1, 21)
    ax2.set_ylim(-1.5, 1.5)
    ax2.set_ylabel('[Fe/H]')

    ax3 = fig.add_subplot(313)
    plt.setp(ax3, xticks=[0, 5, 10, 15, 20], xticklabels=['No cut', 5, 10, 15, 20])
    ax3.errorbar(df.EW, df.vt, yerr=df.vterr, fmt='o', color=c[2])
    ax3.fill_between([0, 20], [1.90-1.08]*2, [1.90+1.08]*2, color='k', alpha=0.2, edgecolor='w')
    ax3.hlines(1.90, 0, 20, linestyles='--')
    ax3.set_xlim(-1, 21)
    ax3.set_ylim(1, 5)
    ax3.set_yticks(range(1, 6))
    ax3.set_xlabel(r'EW cuts [m$\AA$]')
    ax3.set_ylabel(r'$\xi_\mathrm{micro}$ [km/s]')

    plt.tight_layout()
    plt.show()
    # plt.savefig('figures/HD20010_parameters_cuts.pdf')


def fig_synthesis():
    """
    Plot of solar spectrum and a synthesis
    """

    def _plotting(synth, nr=3):
        with open(synth, 'r') as f:
            line = f.readline()
            rows = int(line.split('=')[1].strip(' '))

        data = np.zeros((nr, rows, 2))
        with open(synth, 'r') as f:
            i = -1
            for j, line in enumerate(f):
                if line.startswith('the'):
                    row = 0
                    pass
                elif line.startswith('start'):
                    i += 1
                else:
                    w = float(line[0:12].strip(' '))
                    f = float(line[13::].strip(' '))
                    data[i][row] = [w, f]
                    row += 1
        return data

    data = _plotting('figures/synth.asc', nr=3)
    data1 = data[0][:, 1]

    lines = np.loadtxt('figures/lines80.dat',
                       dtype={'names': ('element', 'w', 'excit', 'gf'),
                              'formats': ('S4', 'f4', 'f4', 'f4')},
                       comments='#', delimiter=',', usecols=(0, 1, 2, 3))

    obs = np.loadtxt('figures/15451.979.asc')
    obs[:, 1] = obs[:, 1]/np.median(obs[:, 1])

    # Get the element from sun.par (line 11 only!!)
    with open('figures/sun.par', 'r') as par:
        for i in range(24):
            par.readline()
        N_elements = int(par.readline().split(' ')[1]) - 1
        par.readline()

    # Setting the format
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.set_xticklabels([])
    ax1.set_ylabel('Normalized flux', fontsize=24)
    ax1.set_xlim(15449, 15455)
    ax1.set_ylim(0.8, 1.05)
    ax2 = fig.add_subplot(212)
    ax2.set_xlim(15449, 15455)
    ax2.xaxis.set_major_formatter(x_formatter)
    ax2.set_xlabel(r'$\lambda$ Angstrom', fontsize=24)
    ax2.set_ylabel('Residuals', fontsize=24)

    # The first plot
    ax1.plot(obs[:, 0], obs[:, 1], '-k', lw=4, alpha=0.6,
             label='Observed spectrum')
    for i, met in zip(range(3), ('0.20', '0.00', '-0.20')):
        lbl = 'Fe abundance: %s' % met
        ax1.plot(data[i][:, 0], data[i][:, 1], label=lbl)
    for line in lines:
        if line[0].startswith('Fe'):
            ax1.vlines(line[1], 0.8, 1.05, alpha=0.3)
    # ax1.legend(frameon=False, loc='best')

    # The second plot
    for i in range(3-1):
        ax2.plot(data[i+1][:, 0], data[i+1][:, 1] - data1)
    ax2.legend((r'$\Delta_{21}$', r'$\Delta_{31}$'), loc='best', frameon=False)
    plt.tight_layout()
    # plt.savefig('figures/synthetic_spectrum.pdf', format='pdf')
    plt.show()


def main():
    """Main function
    """
    # fig_EWvsEP()
    # fig_solarspectrum()
    # fig_solarparams()
    fig_HD20010_parameters()
    # fig_synthesis()


if __name__ == '__main__':
    main()
