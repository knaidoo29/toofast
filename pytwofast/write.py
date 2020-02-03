import numpy as np
import subprocess


def write_paramfile(data_fname, rand_fname, out_fname, param_fname, mode, calc_DD,
                    calc_DR, calc_RR, calc_xi, uselog, useweight, minimum, maximum, numbins, mu_bins):
    """
    Writes a parameter file to be read by twofast.

    Parameters
    ----------
    data_fname : str
        Name of the data file (with full path).
    rand_fname : str
        Name of the randoms file (with full path).
    out_fname : str
        Name of the output file (with full path).
    param_fname : str
        Name of the parameter file to be created (with full path).
    mode : str
        The mode the two point is calculated in 2D, 3D or shell.
    calc_DD : bool
        Auto pair distance for the data.
    calc_DR : bool
        Cross pair distance between data and randoms.
    calc_RR : bool
        Auto pair distance for the randoms.
    calc_xi : bool
        Directly calculate the correlation function using Landy-Szalay estimator.
    uselog : bool
        Where to use logarithmic bins.
    useweight : bool
        Define weight for each point.
    minimum : float
        Minimum distance.
    maximum : float
        Maximum distance.
    numbins : int
        Number of bins to use for auto and cross pair counts.
    mu_bins : int
        Number of bins for multipole expansion of two point.
    """
    subprocess.call('touch '+param_fname, shell=True)
    paramfile = open(param_fname, 'w')
    if data_fname is not None:
        paramfile.write("data_file       " + data_fname + " %\n")
    else:
        paramfile.write("data_file       %\n")
    if rand_fname is not None:
        paramfile.write("rand_file       " + rand_fname + " %\n")
    else:
        paramfile.write("rand_file       %\n")
    paramfile.write("out_file        " + out_fname + " %\n")
    paramfile.write("mode            " + mode + " %\n")
    if calc_DD == True:
        paramfile.write("calc_DD         yes %\n")
    else:
        paramfile.write("calc_DD         no %\n")
    if calc_DR == True:
        paramfile.write("calc_DR         yes %\n")
    else:
        paramfile.write("calc_DR         no %\n")
    if calc_RR == True:
        paramfile.write("calc_RR         yes %\n")
    else:
        paramfile.write("calc_RR         no %\n")
    if calc_xi == True:
        paramfile.write("calc_xi         yes %\n")
    else:
        paramfile.write("calc_xi         no %\n")
    if uselog == True:
        paramfile.write("uselog          yes %\n")
    else:
        paramfile.write("uselog          no %\n")
    if useweight == True:
        paramfile.write("useweight       yes %\n")
    else:
        paramfile.write("useweight       no %\n")
    paramfile.write("minimum         "+str(minimum)+" %\n")
    paramfile.write("maximum         "+str(maximum)+" %\n")
    paramfile.write("numbins         "+str(numbins)+" %\n")
    paramfile.write("mu_bins         "+str(mu_bins)+" %\n")
    paramfile.close()
