# Introduction
The document describes how to calculate Heisenberg exchange parameters for a magnetic compound using Green's functions formalism.

The calculation consists of 3 independent parts:

- Self consistent spin-polarized DFT calculation. Here and below we imply that you use Quantum ESPRESSO (QE) package for such calculation.

- Generation of a model Hamiltonian in Wannier functions basis. The basis should include *d*-states of magnetic ions and sometimes *p*-states of the nearest ligands. This could be done using wannier_ham.x code from QE postprocessing tools. The Hamiltonian should be stored in the [AMULET code](http://amulet-code.org) file format. __Please use the latest stable QE version (5.3.0 and above) for the Hamiltonian production.__

- Calculation of the exchange parameters using exchanges code. The code uses as input only two files: _system.am_ and _hamilt.am_ from the previous step.

# Step 1: Self consistent DFT calculation
Use Quantum ESPRESSO as usual. To obtain a reliable exchange parameters on the last step one should start from a reliable electronic structure on the first step. Perform a spin-polarized calculation on a regular k-points grid within full Brillouin zone. Reciprocal space integration should be performed with gaussian smearing. Use the following keys:
```
nspin = 2
nosym = .true.
noinv = .true.
occupations = 'smearing'
degauss = set a reasonable value for your system
```
One is allowed to perform scf calculation in parallel mode, but `wf_collect = .true.` key must be set, since wannier_ham.x and exchanges.x codes are serial only.

# Step 2: Hamiltonian in Wannier functions basis generation.
The theoretical background for Wannier functions generation and Hamiltonian production procedure is described in [EPJB 65, 91 (2008)](http://www.springerlink.com/index/10.1140/epjb/e2008-00326-3). There is also an example in QE distributive (in `PP/examples/WannierHam_example` dir).

Before you start the model Hamiltonian generation, you should know a symmetry of trial atomic orbitals that will be used for projection (typically these are transition metal d- plus, sometimes, the nearest ligands p-orbitals). And you should know numbers of bands (or energy interval) that you are going to reproduce with the model Hamiltonian. One of the simplest methods to determine the energy interval is to plot the partial densities of states for atomic orbitals with desired symmetry.
__It could be necessary often to include into the basis not only Wannier functions of d-symmetry, but p-orbitals like Wannier functions also. Please note, that exchange interaction parameters are computed on the next step only for d-states. Nevertheless, a superexchange via ligand's p-orbitals is taken into account with nondiagonal elements of Green's functions matrix.__

Structure of the input file for wannier_ham.x (should be sent to stdin):
```
========================================================================
NAMELIST: &INPUTPP
   +--------------------------------------------------------------------
   Variable:       prefix
   
   Type:           CHARACTER
   Default:        ' '
   Description:    as usual
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       outdir
   
   Type:           CHARACTER
   Default:        ' '
   Description:    as usual
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       nwan
   
   Type:           INTEGER
   Default:        ' '
   Description:    Number of Wannier functions
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       use_energy_int
   
   Type:           LOGICAL
   Default:        '.FALSE.'
   Description:    If .true. bands will be defined not by numbers, but by
           energy range (in eV)!
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       plot_bands
   
   Type:           LOGICAL
   Default:        '.FALSE.'
   Description:    If .true. bands structures of original and model hamiltonian
                   will be outputted in gnuplot format for comparison
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       form
   
   Type:           character
   Default:        ' '
   Description:    Hamiltonian file format. Should be set to 'amulet' to obtain
                   the Hamiltonian in the AMULET format suitable for
                   exchanges calculation
   +--------------------------------------------------------------------
========================================================================
CARD: WANNIER_AC

   Definition of trial atomic functions and bands for Wannier generation
   
   /////////////////////////////////////////
   // Syntax:                             //
   /////////////////////////////////////////
   
    Wannier# 1 bands_from bands_to
    atom iatom
    l m
    Wannier# 2 bands_from bands_to
    atom iatom
    l m
    ...
    Spin#2:
    Wannier# 1 bands_from bands_to
    atom iatom
    l m
    Wannier# 2 bands_from bands_to
    atom iatom
    l m
    ...
   
   /////////////////////////////////////////

   ! Please note the Spin#2: block! It should exist (and repeat a Spin#1 block
    to obtain the spin-polarized Hamiltonian).
   
   DESCRIPTION OF ITEMS: 
      +--------------------------------------------------------------------
      Variables:      bands_from, bands_to
      
      Type:           REAL or INTEGER
      Description:    Defines Bloch functions subspace for projection
                      procedure. If use_energy_interval=.true. these are
                      energy values in eV. Otherwise these are bands numbers.
      +--------------------------------------------------------------------
      +--------------------------------------------------------------------
      Variables:      iatom
      
      Type:           INTEGER
      Description:    Number of site on that Wannier function centered
      +--------------------------------------------------------------------
      +--------------------------------------------------------------------
      Variables:      l
      
      Type:           CHARACTER
      Description:    Angular channel for trial wavefunction. 's', 'p' or 'd'
      +--------------------------------------------------------------------
      +--------------------------------------------------------------------
      Variables:      m
      
      Type:           INTEGER
      Description:    Magnetic quantum number of trial orbital (from 1 to 5
                      for d-orbitals, from 1 to 3 for p-orbitals)
      +--------------------------------------------------------------------
===END OF CARD==========================================================
```

Run `wannier_ham.x < your_input_file` on 1 processor in the directory with your scf calculation.

Please check obtained Wannier functions occupations. They should be reliable. If you obtain a lot of `wrong orthogonalization` warnings in output, it means that something is wrong with the energy bands selection.

As a result of Hamiltonian generation procedure at least four files will be produced:
```
hamilt.am
system.am
original_bands.dat
wannier_bands.dat
```
The first two contains the Hamiltonian and some data about crystal structure of the compound __(If you didn't obtain `*.am` files - check `form = 'amulet'` parameter in the input file for wannier_ham.x)__. The last two contains eigenvalues of the full and model Hamiltonians. These files could be plotted with gnuplot and results should coincide within the energy window of interest. The coincidence of the eigenvalues is the main indicator of successful projection procedure.

#Step 3: Exchange parameters calculation with exchanges.x
## Installation
The code requires a Fortran compiler (gfortran is ok) and LAPACK installed in the system. Open `Makefile` with your favorite editor and set `FC` and `LIBS` variables by hands. Than just say `make` inside the directory with source code.

`exchanges.x` file should be generated as a result.

## Usage

The theoretical background for exchanges calculation is described in [Phys. Rev. B 91, 224405 (2015)](http://link.aps.org/doi/10.1103/PhysRevB.91.224405). 

There are few examples of calculations in `examples` dir.

Most of input parameters are set by default and there is no need to change them is most cases. The main thing that should be tuned - `distance` parameter for the nearest neighbors search.

Input file for `exchanges.x` has similar to `wannier_ham.x` Fortran namelists form:
```
========================================================================
NAMELIST: &exchanges
   +--------------------------------------------------------------------
   Variable:       emin
   
   Type:           float
   Default:        -30
   Description:    minimum energy for Green's function integration 
                   relative to Fermi energy in eV
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       emax
   
   Type:           float
   Default:        0.05
   Description:    maximum energy for Green's function integration 
                   relative to Fermi energy in eV
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       height
   
   Type:           float
   Default:        0.5
   Description:    height of the integration contour in eV
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       nz1, nz2, nz3
   
   Type:           integer
   Default:        150,350,150
   Description:    Number of division for integration contour
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       distance
   
   Type:           float
   Default:        0.8
   Description:    Radius of the sphere for the nearest neighbors search
                   in alat
   +--------------------------------------------------------------------
   +--------------------------------------------------------------------
   Variable:       mode
   
   Type:           character
   Default:        'distance'
   Description:    There are two modes of the nearest neihbors search:
                   'distance' and 'list'. In the distance mode the code
                   will calculate exchanges for all d-atoms within the
                   'distance' sphere. In the list mode the code will 
                   consider only atoms (inside the 'distance sphere')
                   listed in ATOMS_LIST card (see below)
   +--------------------------------------------------------------------
========================================================================
========================================================================
CARD: ATOMS_LIST

   Neccessary only if mode='list'
   
   /////////////////////////////////////////
   // Syntax:                             //
   /////////////////////////////////////////
   
    number_of_atoms
    atom_x atom_y atom_z (in crystal coordinates)
    atom_x atom_y atom_z
    atom_x atom_y atom_z
   
   /////////////////////////////////////////
========================================================================
```

# License
The code is distributed under BSD License

# Author
dmitry@korotin.name