#!/usr/bin/env python

import itertools
import math
from datetime import datetime
import errno
import os
import numpy as np


iter_no = 100

mu_d = 1.8550  # D
q_e = 1.0  # Electron charge
convergence_param_d = 1e-4  # D
alpha_prime_ang = 1.45  # Angstrom^3
epsilon_naught = 8.85418782e-12  # F/m
four_pi_epsilon_naught = epsilon_naught * math.pi * 4

# Conversion parameters
debye_conv = 3.33564e-30  # Cm/D
Coulomb_conv = 1.60217646e-19  # C/e
angstrom_conv_cubic = 1e-30  # m^3/Ang^3
angstrom_conv = 1e-10  # m/Ang
electron_volt_conv = 1.602176565e-19  # J/eV


# File names and output directory
outdir = "q3_output"
Field_filename = "Field"
U_ind_out_filename = "U_ind"
mu_ind_out_filename = "mu_ind"
dipoles_out_filename = "dipoles"
del_mu_ind_out_filename = "delta_mu_ind"
del_mu_ind_check_out_filename = "delta_mu_ind_check"


# Timer
startTime = datetime.now()

# Change units from Debye to Coulomb-Meter
mu = mu_d * debye_conv
convergence_param = convergence_param_d * debye_conv

# Change units from electron charge to Coulombs
q = q_e * Coulomb_conv

# Change units from Angstrom^3 to Meter^3
alpha_prime = alpha_prime_ang * angstrom_conv_cubic


# Create output directory function
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# Make the output directory
mkdir_p(outdir)


# Create filename function
def filename_change(path, filename, distance):

    temp = '/'.join([path, filename])
    # d = str(distance)
    d = "{0:.1f}".format(distance)
    temp2 = '_'.join([temp, d])
    current_filename = '.'.join([temp2, "txt"])

    return current_filename

print "\nOutput\nConvergence requirement (Cm):", convergence_param, "\n"

site_no = [1, 2, 3, 4, 5]

radiusList = [2.5, 3.0, 3.5]
# radiusList = [2.5]

for rad in range(len(radiusList)):

    # Generate names for output files
    Fieldout = filename_change(outdir, Field_filename, radiusList[rad])
    U_ind_out = filename_change(outdir, U_ind_out_filename, radiusList[rad])
    mu_ind_out = filename_change(outdir, mu_ind_out_filename, radiusList[rad])
    dipoles_out = filename_change(outdir, dipoles_out_filename, radiusList[rad])
    del_mu_ind_out = filename_change(outdir, del_mu_ind_out_filename, radiusList[rad])
    del_mu_ind_check_out = filename_change(outdir, del_mu_ind_check_out_filename, radiusList[rad])

    print "radius (Angstroms)", radiusList[rad]

    # Change units from Angstrom to Meter
    radius = radiusList[rad] * angstrom_conv

    qpos = -radius

    position = []

    # setup dipole positions
    for i in range(len(site_no)):
        position.append(radius * i)

    # Intiate variables
    U_ind = [0, 0, 0, 0, 0]
    mu_ind = [0, 0, 0, 0, 0]
    Field = [0, 0, 0, 0, 0]
    del_mu_ind = [1, 1, 1, 1, 1]
    U_elect = [0, 0, 0, 0, 0]
    del_mu_ind_check = []
    dipoles = []

    ion_only_mu_ind = [0, 0, 0, 0, 0]
    ion_only_dipoles = []
    ion_only_Field = [0, 0, 0, 0, 0]

    for i in itertools.repeat(mu, 5):
        dipoles.append(i)
        ion_only_dipoles.append(i)

    for i in itertools.repeat("No", 5):
        del_mu_ind_check.append(i)


    # Calculate the permanent electrostatic energy and write out to screen
    for dipole in range(len(dipoles)):

        # charge-dipole interactions
        U_elect[dipole] = (- q * mu) / (four_pi_epsilon_naught * (position[dipole] - qpos) ** 2)

        # dipole-dipole interactions
        for partner in range(len(dipoles)):

            if partner != dipole:  # Count all interaction pairs but no self-interaction
                dist = abs(position[dipole] - position[partner])

                U_elect[dipole] += - (2 * mu ** 2) / ((dist ** 3) * four_pi_epsilon_naught)


    # Calculate the ion only induced dipole Field and induced dipole
    for dipole in range(len(dipoles)):
        ion_only_Field[dipole] = q / (four_pi_epsilon_naught * ((position[dipole] - qpos) ** 2))

        ion_only_mu_ind[dipole] = alpha_prime * four_pi_epsilon_naught * ion_only_Field[dipole]

        ion_only_dipoles[dipole] = mu + ion_only_mu_ind[dipole]

    # Iterative induced dipole loop (using with to open output files)
    with open(Fieldout, 'w') as Fieldfile, \
            open(U_ind_out, 'w') as U_ind_file, \
            open(mu_ind_out, 'w') as mu_ind_file, \
            open(dipoles_out, 'w') as dipoles_file, \
            open(del_mu_ind_out, 'w') as del_mu_ind_file, \
            open(del_mu_ind_check_out, 'w') as del_mu_ind_check_file:


        for num in range(0, iter_no):

            for dipole in range(len(dipoles)):  # Calculate the current field at each dipole

                # Field from the ion
                Field[dipole] = q / (four_pi_epsilon_naught * ((position[dipole] - qpos) ** 2))

                # Field from the dipoles
                for partner in range(len(dipoles)):

                    if partner != dipole:  # Count all interaction pairs but no self-interaction
                        dist = abs(position[dipole] - position[partner])
                        Field[dipole] += (2 * dipoles[partner]) / ((dist ** 3) * four_pi_epsilon_naught)


                # Calculate the induced dipole, its energy and it's change per cycle
                U_ind[dipole] = -0.5 * (Field[dipole]**2) * alpha_prime * four_pi_epsilon_naught

                del_mu_ind[dipole] = -mu_ind[dipole]

                mu_ind[dipole] = alpha_prime * four_pi_epsilon_naught * Field[dipole]

                dipoles[dipole] = mu + mu_ind[dipole]

                del_mu_ind[dipole] += mu_ind[dipole]

                if del_mu_ind[dipole] > convergence_param:
                    del_mu_ind_check[dipole] = "no"
                else:
                    del_mu_ind_check[dipole] = "Converged"


            #  Write out calculated numbers
            Fieldfile.write("%s\n" % Field)
            U_ind_file.write("%s\n" % U_ind)
            mu_ind_file.write("%s\n" % mu_ind)
            dipoles_file.write("%s\n" % dipoles)
            del_mu_ind_file.write("%s\n" % del_mu_ind)
            del_mu_ind_check_file.write("%s\n" % del_mu_ind_check)


    #  Write out values
    print "\nSI units:"
    print "U_elect components (J)", U_elect
    print "U_elect (J)", sum(U_elect)
    print "U_ind components (J)", U_ind
    print "U_ind (J)", sum(U_ind)
    print "U_total (J)", sum(U_ind) + sum(U_elect)
    print "mu_ind (Cm)", mu_ind
    print "total dipoles (Cm)", dipoles
    print "ion only induced dipoles (Cm)", ion_only_mu_ind
    print "ion only total dipoles (Cm)", ion_only_dipoles

    U_elect_ev = [0, 0, 0, 0, 0]
    U_ind_ev = [0, 0, 0, 0, 0]
    mu_ind_debye = [0, 0, 0, 0, 0]
    dipoles_debye = [0, 0, 0, 0, 0]
    ion_only_mu_ind_debye = [0, 0, 0, 0, 0]
    ion_only_dipoles_debye = [0, 0, 0, 0, 0]


    # Convert into Debye and eV
    for dipole in range(len(dipoles)):
        U_elect_ev[dipole] = U_elect[dipole] / electron_volt_conv
        U_ind_ev[dipole] = U_ind[dipole] / electron_volt_conv
        mu_ind_debye[dipole] = mu_ind[dipole] / debye_conv
        dipoles_debye[dipole] = dipoles[dipole] / debye_conv
        ion_only_mu_ind_debye[dipole] = ion_only_mu_ind[dipole] / debye_conv
        ion_only_dipoles_debye[dipole] = ion_only_dipoles[dipole] / debye_conv


    U_elect_ev_out = np.array(U_elect_ev)
    U_ind_ev_out = np.array(U_ind_ev)
    mu_ind_debye_out = np.array(mu_ind_debye)
    dipoles_debye_out = np.array(dipoles_debye)
    ion_only_mu_ind_debye_out = np.array(ion_only_mu_ind_debye)
    ion_only_dipoles_debye_out = np.array(ion_only_dipoles_debye)
    np.set_printoptions(precision=5)

    print "\nnon-SI units:"
    print "U_elect components (eV)", U_elect_ev_out
    print "U_elect (eV)", np.sum(U_elect_ev_out)
    print "U_ind components (eV)", U_ind_ev_out
    print "U_ind (eV)", np.sum(U_ind_ev_out)
    print "U_total (eV)", np.sum(U_ind_ev_out) + np.sum(U_elect_ev_out)
    print "mu_ind (Debye)", mu_ind_debye_out
    print "total dipoles (Debye)", dipoles_debye_out
    print "ion only induced dipoles (Debye)", ion_only_mu_ind_debye_out
    print "ion only total dipoles (Debye)", ion_only_dipoles_debye_out, "\n\n"

print "\nExecution time", datetime.now() - startTime
