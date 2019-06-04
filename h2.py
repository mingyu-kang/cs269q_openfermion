from openfermion.hamiltonians import MolecularData, load_molecular_hamiltonian
from openfermion.ops import FermionOperator
from openfermion.transforms import bravyi_kitaev
from openfermion.utils import hermitian_conjugated
from openfermion.ops import QubitOperator
from forestopenfermion import pyquilpauli_to_qubitop, qubitop_to_pyquilpauli
from pyquil.paulis import PauliSum
from pyquil.api import WavefunctionSimulator
from scipy.optimize import minimize
from pyquil import Program
from pyquil.gates import *

import numpy as np
import functools
import matplotlib.pyplot as plt

sim = WavefunctionSimulator()

numQubit = 4

class Track:
    def __init__(self):
        self.cnt = 0

def solve_vqe(hamiltonian: PauliSum, numLayer):
    # Construct a variational quantum eigensolver solution to find the lowest
    # eigenvalue of the given hamiltonian
    theta_init = np.random.rand(2*numQubit*numLayer) * 2 * np.pi

    track = Track()

    def inc(t, xk):
        t.cnt += 1

    def ansatz_energy(theta_vec):
        p = None
        p = Program()
        for j in range(numLayer):
            for i in range(numQubit):
                p += RX(theta_vec[2*numQubit * j + 2 * i], i)
                p += RZ(theta_vec[2*numQubit * j + 2 * i + 1], i)
            for i in range(numQubit - 1):
                p += CNOT(i, i + 1)

        energy = sim.expectation(p, hamiltonian).real
        return energy

    theta_answer = minimize(ansatz_energy, theta_init, method='L-BFGS-B', callback=functools.partial(inc, track)).x
    return ansatz_energy(theta_answer), track.cnt


def get_ground_energy(interaction_hamil, numLayer):
    fermionop_hamil = FermionOperator()
    for key in interaction_hamil:
        value = interaction_hamil[key]
        fermionop_hamil += FermionOperator(term=key, coefficient=value)

    qubitop_hamil = bravyi_kitaev(fermionop_hamil)
    pauliop_hamil = qubitop_to_pyquilpauli(qubitop_hamil)

    return solve_vqe(pauliop_hamil, numLayer)

def e_layers(bond_lengths):
    basis = 'sto-3g'
    multiplicity = 1  # 2S+1

    data = []

    for layer in [1,2,3,4,5]:
        gelist = []
        for bond_length in bond_lengths:
            geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
            description = str(round(bond_length, 2))
            h2_interaction_hamil = load_molecular_hamiltonian(geometry,
                basis,
                multiplicity,
                description,
                n_active_electrons=None,
                n_active_orbitals=None)

            ge = get_ground_energy(h2_interaction_hamil, layer)
            gelist.append(ge)
            print('bond length: ', round(bond_length, 2), ' ground state energy: ', ge)
        data.append(gelist)

    return data

def plot_data(bond_lengths, data):
    basis = 'sto-3g'
    multiplicity = 1
    bond_length_interval = 0.1
    n_points = 25

    # Generate molecule at different bond lengths.
    hf_energies = []
    fci_energies = []

    for bond_length in bond_lengths:
        description = str(round(bond_length, 2))
        #    print(description)
        geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
        molecule = MolecularData(
            geometry, basis, multiplicity, description=description)

        # Load data.
        molecule.load()

        hf_energies += [molecule.hf_energy]
        fci_energies += [molecule.fci_energy]

    plt.figure(0)
    plt.plot(bond_lengths, fci_energies, 'x-')
    plt.plot(bond_lengths, [e for e, _ in data], 'o-')
    plt.ylabel('Energy in Hartree')
    plt.xlabel('Bond length in angstrom')
    plt.show()

if __name__ == "__main__":
    bond_lengths = np.linspace(0.3, 2.5, 23)
    data = e_layers(bond_lengths)
    plot_data(data)
