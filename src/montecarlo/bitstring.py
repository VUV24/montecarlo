import numpy as np
import math      
import copy as cp       


class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __repr__(self):
        out = ""
        for i in self.config:
            out += str(i)
        return out

    def __eq__(self, other):        
        return all(self.config == other.config)
    
    def __len__(self):
        return len(self.config)

    def on(self):
        """
        Return number of bits that are on
        """
        return np.count_nonzero(self.config)

    def off(self):
        """
        Return number of bits that are off
        """
        return self.N - self.on()

    def flip_site(self,i):
        """
        Flip the bit at site i
        """
        self.config[i] = 1 - self.config[i]
    
    def integer(self):
        """
        Return the decimal integer corresponding to BitString
        """
        dec = 0
        for i in range(len(self.config)):
            if self.config[i] == 1:
                dec += pow(2, self.N - 1 - i)

        return dec

 

    def set_config(self, s:list[int]):
        """
        Set the config from a list of integers
        """
        self.config = np.array(s, dtype = int)

    def set_integer_config(self, dec:int):
        """
        convert a decimal integer to binary
    
        Parameters
        ----------
        dec    : int
            input integer
            
        Returns
        -------
        Bitconfig
        """
        for i in range(len(self.config)):
            self.config[self.N - 1 - i] = dec % 2
            dec //= 2