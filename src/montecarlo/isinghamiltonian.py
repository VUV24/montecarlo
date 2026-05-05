import numpy as np
import networkx as nx
import random
import scipy
from .bitstring import BitString

class IsingHamiltonian:
    """
    Represents an Ising Hamiltonian on a graph.
    """
    def __init__(self, G: nx.Graph):
        self.G = G
        self.mus = np.zeros(G.number_of_nodes())
        
    def energy(self, bs: BitString):
        """
        Compute the energy of a given spin configuration.
        """
        G = self.G
        
        energy = 0.0
        right = 0.0
        for i in G.edges:
            si = (bs.config[i[0]] * 2) - 1
            sj = (bs.config[i[1]] * 2) - 1
            energy += G.edges[i]["weight"] * si * sj
            
        for i in range(len(self.mus)):
            right += self.mus[i] * ((bs.config[i] * 2) - 1)

        energy += right
        return energy

    def compute_average_values(self, T: float):
        """
        Compute thermodynamic averages by summing over all configurations.
        """
        E  = 0.0
        M  = 0.0
        Z  = 0.0
        EE = 0.0
        MM = 0.0

        bs = BitString(self.G.number_of_edges())

        for i in range(2**bs.N):
            bs.set_integer_config(i)
            exp = -1 * (1/T) * self.energy(bs)
            Z += np.e ** exp
        
        for i in range(2**bs.N):
            bs.set_integer_config(i)
            exp = -1 * (1/T) * self.energy(bs)
            prob = (np.e**exp) / Z
            E += self.energy(bs) * prob
            numOn = bs.on()
            numOff = bs.off()
            m_alpha = numOn - numOff
            M += m_alpha * prob
            EE += (self.energy(bs)**2) * prob
            MM += (m_alpha**2) * prob

        HC = (EE - (E**2)) / (T**2)
        MS = (MM - (M**2)) / T
        
        return E, M, HC, MS 
    
    def set_mu(self, mus: np.array):
        """
        Set the local mus for each site.
        """
        self.mus = mus
