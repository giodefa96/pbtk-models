from abc import ABC, abstractmethod

class Molecule(ABC):
    @abstractmethod
    def __init__(self, V_max, K_m, P_f, P_l, P_s, P_r, P_b):
        # Metabolism parameters
        self.V_max = V_max  # Maximal velocity of metabolism
        self.K_m = K_m  # Michaelis–Menten affinity

        # Molecule interaction with tissues
        self.P_f = P_f  # Fat:blood partition coefficient
        self.P_l = P_l  # Liver:blood partition coefficient
        self.P_s = P_s  # Poorly perfused tissue:blood partition coefficient
        self.P_r = P_r  # Richly perfused tissue:blood partition coefficient
        self.P_b = P_b  # Blood:air partition coefficient
