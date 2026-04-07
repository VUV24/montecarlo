from .isinghamiltonian import IsingHamiltonian

class MonteCarlo:
    def __init__(self, ham: IsingHamiltonian):
        self.ham = ham

    def run(self, T: int, n_samples: int, n_burn: int):
        G = self.ham.G
        for i in range(G.number_of_nodes()):
            G = self.ham.G