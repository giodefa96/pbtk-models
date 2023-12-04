from .molecules import Molecule

class Perchloroethylene(Molecule):
    def __init__(self, V_max=1.66, K_m=0.55, P_f=56.72, P_l=4.64, P_s=1.54, P_r=4.64, P_b=18):
        super().__init__(V_max, K_m, P_f, P_l, P_s, P_r, P_b)
