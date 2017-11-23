#
# @BEGIN LICENSE
#
# syshfwplug by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util
import numpy as np
import sys,struct
import scipy
from scipy import optimize


def run_syshfwplug(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    syshfwplug can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('syshfwplug')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.core.set_local_option('MYPLUGIN', 'PRINT', 1)

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    syshfwplug_wfn = psi4.core.plugin('syshfwplug.so', ref_wfn)

    return syshfwplug_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['syshfwplug'] = run_syshfwplug




def localize_PM(W):
    C     = np.asarray(W.Ca_subset("AO","OCC"))
    def transform_and_calc(para):
        X = np.ones((int(len(para)/2),int(len(para)/2)))
        x = np.triu_indices(int(len(para)/2),1)
        X[x] = para
        X+=-X.transpose()
        Uni = scipy.linalg.expm(X)
        L = np.matmul(C,Uni)
        PM_spread = calculate_PM(W,L)
        #print np.sum(PA**2)
        return 1/PM_spread
        
    U = np.random.randint(0,20,int((C.shape[1]-1)*C.shape[1]/2))/100
    results = scipy.optimize.minimize(transform_and_calc,U,method='BFGS')
    X = np.ones((int(len(U)/2),int(len(U)/2)))
    x = np.triu_indices(int(len(U)/2),1)
    X[x] = results.x
    X+=-X.transpose()
    Uni = scipy.linalg.expm(X)
    L = np.matmul(C,Uni)
    return results,L





def localize_boys(W):
    C     = np.asarray(W.Ca_subset("AO","OCC"))
    def transform_and_calc(para):
        X = np.ones((int(len(para)/2),int(len(para)/2)))
        x = np.triu_indices(int(len(para)/2),1)
        X[x] = para
        X+=-X.transpose()
        Uni = scipy.linalg.expm(X)

        L = np.matmul(C,Uni)
        boys_spread = calculate_boys(W,L)
        return 1/boys_spread

    U = np.zeros(int((C.shape[1]-1)*C.shape[1]/2))/100
    results = scipy.optimize.minimize(transform_and_calc,U,method='BFGS')

    #construct the final antisymmetric matrix
    X = np.ones((int(len(U)/2),int(len(U)/2)))
    x = np.triu_indices(int(len(U)/2),1)
    X[x] = results.x
    X+=-X.transpose()

    Uni = scipy.linalg.expm(X)
    L = np.matmul(C,Uni)
    return results,L



def calculate_PM(W,L):
    basis = W.basisset()  # we need the basis
    S     = np.asarray(W.S()) # Overlap
    natom = W.molecule().natom()
    basm = np.array([basis.function_to_center(x) for x in range(basis.nbf())])
    LS = np.matmul(S,L)
    PA    = np.zeros((natom,L.shape[1]))
    for i in range(L.shape[1]):
        for j in list(set(basm)):
            for k in [c for c,x in enumerate(basm) if x==j]: # k contains the indices,locates at center i
                PA[j][i] += L[k][i]*LS[k][i]
    return np.sum(PA**2)



def calculate_boys(W,L):
    mints = psi4.core.MintsHelper(W.basisset())
    Dx = np.asarray(mints.ao_dipole()[0])
    Dy = np.asarray(mints.ao_dipole()[1])
    Dz = np.asarray(mints.ao_dipole()[2])

    DXmo = np.matmul(np.matmul(L.transpose(),Dx), L)
    DYmo = np.matmul(np.matmul(L.transpose(),Dy), L) 
    DZmo = np.matmul(np.matmul(L.transpose(),Dz), L)
    
    return (np.sum(np.diag(DXmo)**2)+np.sum(np.diag(DYmo)**2) + np.sum(np.diag(DZmo)**2))
    

    

def write_wavefunction_file(filename,eps,Mos):
    f = open(filename,'wb')
    nroao = int(Mos.shape[0])
    f.write(struct.pack('i', nroao ))
    f.write(struct.pack(str(nroao)+'d',*eps))
    f.write(struct.pack(str(nroao**2)+'d',*Mos.transpose().reshape(nroao**2)))
    f.close()




  


def exampleFN(W):
    """
    resB,LB = localize_boys(W)
    C= np.asarray(W.Ca_subset("AO","OCC"))
    #print calculate_boys(W,C)
    #print calculate_boys(W,LB)
    loc = psi4.core.Localizer.build("BOYS",W.basisset(),W.Ca_subset("AO","OCC"))
    loc.localize()

    #print calculate_boys(W,np.asarray(loc.L))
    loc = psi4.core.Localizer.build("PIPEK_MEZEY",W.basisset(),W.Ca_subset("AO","OCC"))
    loc.localize()
    #print calculate_PM(W,C)
    resPM,LPM = localize_PM(W)
    #print calculate_PM(W,LPM)
    #print calculate_PM(W,np.asarray(loc.L))


    C = np.zeros((LPM.shape[0],LPM.shape[0]))
    C[:,0:LPM.shape[1]] = LPM
    CVir = np.asarray(W.Ca_subset("AO","VIR"))
    C[:,LPM.shape[1]:] = CVir
    
    write_wavefunction_file('MyPM.ahfw',np.zeros(C.shape[0]),C)

   
    molden = psi4.core.MoldenWriter(W)
    molden.write(str(u'loc.molden'),psi4.core.Matrix.from_array(C),psi4.core.Matrix.from_array(C),W.epsilon_a(),W.epsilon_a(),W.epsilon_a(),W.epsilon_a(),True)

    print W.basisset().nbf()    
  
    base_wfn = psi4.core.Wavefunction.build(mol , psi4.core.get_global_option('BASIS'))
    basis   = psi4.core.BasisSet.build(mol)
    ribasis = psi4.core.BasisSet.build(mol, 'DF_BASIS_MP2', '', 'RIFIT', psi4.core.get_global_option('BASIS'))

    
    mints = psi4.core.MintsHelper(base_wfn.basisset())

    T = np.asarray(mints.ao_kinetic())
    V = np.asarray(mints.ao_potential())
    ECP = np.asarray(mints.ao_ecp())
    
    H = T+V+ECP
    """
    aux = psi4.core.BasisSet.build(mol, "DF_BASIS_MP2", "", "RIFIT", "aug-cc-pvdz")


    return 0





