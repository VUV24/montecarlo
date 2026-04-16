from .isinghamiltonian import IsingHamiltonian
from .bitstring import BitString
import random
import numpy as np

class MonteCarlo:
    def __init__(self, ham: IsingHamiltonian):
        self.ham = ham

    # n_samples = total num of steps in montecarlo, n_burn = number of configs to skip before recording
    def run(self, T: int, n_samples: int, n_burn: int):

        energies = []
        mags = []

        bs = BitString(self.ham.G.number_of_edges())

        for i in range(n_samples):
            bs.set_integer_config(i)

            for j in range(bs.__len__()):
                if (n_burn < 0):
                    preEner = self.ham.energy(bs)
                    bs.flip_site(j)
                    postEner = self.ham.energy(bs)
                    energyDiff = postEner - preEner

                    if (energyDiff <= 0):
                        energies.append(postEner)
                        mag = bs.on() - bs.off() 
                        mags.append(mag)

                    elif (energyDiff > 0):
                        r = random.random()
                        num = (np.e ** (-1.0 * energyDiff))   
                        if (num > r):
                            mag = bs.on() - bs.off()
                            energies.append(postEner)
                            mags.append(mag)
            n_burn -= 1

        return energies, mags

                
            