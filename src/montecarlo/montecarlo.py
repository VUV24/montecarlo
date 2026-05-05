from .isinghamiltonian import IsingHamiltonian
from .bitstring import BitString
import random
import numpy as np

class MonteCarlo:
    """
    Performs Monte Carlo simulation of an Ising Hamiltonian using Metropolis sampling.
    """
    def __init__(self, ham: IsingHamiltonian):
        self.ham = ham

    # n_samples = total num of steps in montecarlo, n_burn = number of configs to skip before recording
    def run(self, T: int, n_samples: int, n_burn: int):
        """
        Perform the Metropolis algorithm.
        """
        energies = []
        mags = []

        bs = BitString(self.ham.G.number_of_nodes())
        bs.set_integer_config(random.randint(0, (2**bs.__len__()) - 1))

        for i in range(n_samples+n_burn):
            for j in range(bs.__len__()):
                preEner = self.ham.energy(bs)
                bs.flip_site(j)
                postEner = self.ham.energy(bs)
                energyDiff = postEner - preEner

                if (energyDiff <= 0):
                    pass
                elif (energyDiff > 0):
                    r = random.random()
                    num = (np.e ** (-1.0 * energyDiff / T))
                    if (num > r):
                        pass
                    else:
                        bs.flip_site(j)

            if (i >= n_burn):
                energies.append(self.ham.energy(bs))
                mag = bs.on() - bs.off() 
                mags.append(mag)

        return energies, mags