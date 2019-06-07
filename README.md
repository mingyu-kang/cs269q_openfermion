# CS269Q Final Project

This is the repository for CS269Q final project - **Finding Ground State Energy of Molecules with Variational Quantum Eigensolver** by Mingyu Kang and Jaebum Lee

### Requirements

You need to have up-to-date versions of `openfermion`, `forestopenfermion`

Also, you may need to install [psi4 using the guideline here](http://www.psicode.org/psi4manual/1.2/external.html) to check the numerically computed ground state energies. The default version of the openfermion only contains the data of H2. You must run psi4 to generate the data for LiH and BeH2. If you don't want to run the program by yourself, we have also provided pre-computed data for LiH and BeH2 in `data` folder of this repository. You **must** copy the data into openfermion's package directory.

### Files

Those two computes the ground state energy of H2 using hardware ansatz and UCC respectively.

- H2 Hardware Ansatz.ipynb
- H2 UCC.ipynb

Those two computes the ground state energy of LiH using hardware ansatz and UCC respectively.

- LiH Hardware Ansatz.ipynb
- LiH_UCC.ipynb

Those two computes the ground state energy of BeH2 using UCC.

-  BeH2_UCC.ipynb
