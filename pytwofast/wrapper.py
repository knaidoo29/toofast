import numpy as np
import subprocess
import os
from . import write
from . import source

class TwoPoint:

    def __init__(self):
        self.mode_data = None
        self.mode_rand = None
        self.path = os.getcwd()+'/'
        self.minimum = None
        self.maximum = None
        self.numbins = None
        self.mu_bins = None
        self.calc_DD = None
        self.calc_DR = None
        self.calc_RR = None
        self.calc_xi = None
        self.uselog = None
        self.useweight = None
        self.identifier = None

    def setup(self, minimum, maximum, numbins, mu_bins=10, calc_DD=True, calc_DR=True,
              calc_RR=True, calc_xi=True, uselog=False, useweight=False, identifier=None):
        self.minimum = minimum
        self.maximum = maximum
        self.numbins = numbins
        self.mu_bins = mu_bins
        self.calc_DD = calc_DD
        self.calc_DR = calc_DR
        self.calc_RR = calc_RR
        self.calc_xi = calc_xi
        self.uselog = uselog
        self.useweight = useweight
        self.identifier = identifier

    def add_data(self, x=None, y=None, z=None, phi=None, theta=None, weight=None):
        if self.useweight is True and weight is None:
            print("useweight is true but weights have not been assigned, assigning weights to 1.")
            if x is not None:
                weight = np.ones(len(x))
            elif phi is not None:
                weight = np.ones(len(phi))
        if x is not None and y is not None and z is None and phi is None and theta is None:
            self.mode_data = '2D'
            if self.identifier is None:
                if self.useweight is False:
                    np.savetxt(self.path + "temp_data.txt", zip(x, y))
                else:
                    np.savetxt(self.path + "temp_data.txt", zip(x, y, weight))
            else:
                if self.useweight is False:
                    np.savetxt(self.path + "temp_data_"+str(self.identifier)+".txt", zip(x, y))
                else:
                    np.savetxt(self.path + "temp_data_"+str(self.identifier)+".txt", zip(x, y, weight))
        elif x is not None and y is not None and z is not None and phi is None and theta is None:
            self.mode_data = '3D'
            if self.identifier is None:
                if self.useweight is False:
                    np.savetxt(self.path + "temp_data.txt", zip(x, y, z))
                else:
                    np.savetxt(self.path + "temp_data.txt", zip(x, y, z, weight))
            else:
                if self.useweight is False:
                    np.savetxt(self.path + "temp_data_"+str(self.identifier)+".txt", zip(x, y, z))
                else:
                    np.savetxt(self.path + "temp_data_"+str(self.identifier)+".txt", zip(x, y, z, weight))
        elif x is None and y is None and z is None and phi is not None and theta is not None:
            self.mode_data = 'tomo'
            condition = np.where((phi < 0.) | (phi > 2.*np.pi) | (theta < 0.) | (theta > np.pi))[0]
            if len(condition) == 0:
                if self.identifier is None:
                    if self.useweight is False:
                        np.savetxt(self.path + "temp_data.txt", zip(phi, theta))
                    else:
                        np.savetxt(self.path + "temp_data.txt", zip(phi, theta, weight))
                else:
                    if self.useweight is False:
                        np.savetxt(self.path + "temp_data_"+str(self.identifier)+".txt", zip(phi, theta))
                    else:
                        np.savetxt(self.path + "temp_data_"+str(self.identifier)+".txt", zip(phi, theta, weight))
            else:
                print("!!Tomographic Range Error!!")
                print("phi must be given in radians between 0 and 2*PI")
                print("current range for phi = ", phi.min(), phi.max())
                print("theta must be given in radians between 0 and PI")
                print("current range for theta = ", theta.min(), theta.max())
        elif x is None and y is None and z is None and phi is None and theta is None:
            self.mode_data = None
        else:
            self.mode_data = None

    def add_rand(self, x=None, y=None, z=None, phi=None, theta=None, weight=None):
        if self.useweight is True and weight is None:
            print("useweight is true but weights have not been assigned, assigning weights to 1.")
            if x is not None:
                weight = np.ones(len(x))
            elif phi is not None:
                weight = np.ones(len(phi))
        if x is not None and y is not None and z is None and phi is None and theta is None:
            self.mode_rand = '2D'
            if self.identifier is None:
                if self.useweight is False:
                    np.savetxt(self.path + "temp_rand.txt", zip(x, y))
                else:
                    np.savetxt(self.path + "temp_rand.txt", zip(x, y, weight))
            else:
                if self.useweight is False:
                    np.savetxt(self.path + "temp_rand_"+str(self.identifier)+".txt", zip(x, y))
                else:
                    np.savetxt(self.path + "temp_rand_"+str(self.identifier)+".txt", zip(x, y, weight))
        elif x is not None and y is not None and z is not None and phi is None and theta is None:
            self.mode_rand = '3D'
            if self.identifier is None:
                if self.useweight is False:
                    np.savetxt(self.path + "temp_rand.txt", zip(x, y, z))
                else:
                    np.savetxt(self.path + "temp_rand.txt", zip(x, y, z, weight))
            else:
                if self.useweight is False:
                    np.savetxt(self.path + "temp_rand_"+str(self.identifier)+".txt", zip(x, y, z))
                else:
                    np.savetxt(self.path + "temp_rand_"+str(self.identifier)+".txt", zip(x, y, z, weight))
        elif x is None and y is None and z is None and phi is not None and theta is not None:
            self.mode_rand = 'tomo'
            condition = np.where((phi < 0.) | (phi > 2.*np.pi) | (theta < 0.) | (theta > np.pi))[0]
            if len(condition) == 0:
                if self.identifier is None:
                    if self.useweight is False:
                        np.savetxt(self.path + "temp_rand.txt", zip(phi, theta))
                    else:
                        np.savetxt(self.path + "temp_rand.txt", zip(phi, theta, weight))
                else:
                    if self.useweight is False:
                        np.savetxt(self.path + "temp_rand_"+str(self.identifier)+".txt", zip(phi, theta))
                    else:
                        np.savetxt(self.path + "temp_rand_"+str(self.identifier)+".txt", zip(phi, theta, weight))
            else:
                print("!!Tomographic Range Error!!")
                print("phi must be given in radians between 0 and 2*PI")
                print("current range for phi = ", phi.min(), phi.max())
                print("theta must be given in radians between 0 and PI")
                print("current range for theta = ", theta.min(), theta.max())
        elif x is None and y is None and z is None and phi is None and theta is None:
            self.mode_rand = None
        else:
            self.mode_rand = None

    def run(self, output_fname, location='home', processors=4, outmultipole=False):
        if outmultipole is True:
            self.mode = "poly"
        if self.mode_data is not None:
            if self.identifier is None:
                data_fname = self.path + "temp_data.txt"
            else:
                data_fname = self.path + "temp_data_"+str(self.identifier)+".txt"
        else:
            data_fname = None
        if self.mode_rand is not None and self.mode_data == self.mode_rand:
            if self.identifier is None:
                rand_fname = self.path + "temp_rand.txt"
            else:
                rand_fname = self.path + "temp_rand_"+str(self.identifier)+".txt"
        else:
            rand_fname= None
        out_fname = self.path + output_fname
        if self.identifier is None:
            param_fname = self.path + "temp_paramfile.ini"
        else:
            param_fname = self.path + "temp_paramfile_"+str(self.identifier)+".ini"
        mode = self.mode_data
        write.write_paramfile(data_fname, rand_fname, out_fname, param_fname, mode,
                              self.calc_DD, self.calc_DR, self.calc_RR, self.calc_xi,
                              self.uselog, self.useweight, self.minimum, self.maximum, self.numbins)
        mpirun, twopoint, twopoint_mpi = source.get_src(location)
        if processors == 1:
            subprocess.call(twopoint + " " +param_fname, shell=True)
        else:
            if location == 'splinter':
                subprocess.call("mpirun " + twopoint_mpi + " " + param_fname, shell=True)
            else:
                subprocess.call(mpirun + " -n " + str(processors) + " " + twopoint_mpi + " " +param_fname, shell=True)

    def clean(self, remove_temp_files=True):
        if remove_temp_files == True:
            if self.identifier is None:
                subprocess.call("rm " + self.path + "temp_data.txt", shell=True)
                subprocess.call("rm " + self.path + "temp_rand.txt", shell=True)
                subprocess.call("rm " + self.path + "temp_paramfile.ini", shell=True)
            else:
                subprocess.call("rm " + self.path + "temp_data_"+str(self.identifier)+".txt", shell=True)
                subprocess.call("rm " + self.path + "temp_rand_"+str(self.identifier)+".txt", shell=True)
                subprocess.call("rm " + self.path + "temp_paramfile_"+str(self.identifier)+".ini", shell=True)
        self.__init__()
