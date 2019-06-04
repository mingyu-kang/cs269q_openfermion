import os

from openfermion.hamiltonians import MolecularData

from openfermionpsi4 import run_psi4

if __name__ == '__main__':

    # Set chemical parameters.
    basis = 'sto-3g'
    charge = 0
    multiplicity = 1

    # Single point at equilibrium for testing
    spacings = [1.7698]

    # Add points for a full dissociation curve from 0.1 to 3.0 angstroms
    spacings += [0.1 * r for r in range(1, 30)]

    # Set run options
    run_scf = 1
    run_mp2 = 1
    run_cisd = 1
    run_ccsd = 1
    run_fci = 1
    verbose = 1
    tolerate_error = 1

    # Run Diatomic Curve
    for spacing in spacings:
        description = "{:.2f}".format(spacing)
        geometry = [['Be', [0, 0, 0]],
                    ['H', [0, 0, spacing]],
                    ['H', [0, 0, -spacing]]]
        molecule = MolecularData(geometry,
                                 basis,
                                 multiplicity,
                                 charge,
                                 description)

        molecule = run_psi4(molecule,
                            run_scf=run_scf,
                            run_mp2=run_mp2,
                            run_cisd=run_cisd,
                            run_ccsd=run_ccsd,
                            run_fci=run_fci,
                            verbose=verbose,
                            tolerate_error=tolerate_error)
        molecule.save()