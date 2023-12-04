from .molecules import Molecule


class Tetrachloroethane(Molecule):
    """
    1,1,2,2-Tetrachloroethane
    """
    def __init__(self, V_max=0, K_m=0.000001, P_l=2.55, P_f=44.9, P_r=2.55, P_s=2.46, P_k=3.07, P_b=5.17):
        super().__init__(V_max, K_m, P_f, P_l, P_s, P_r, P_b)
        self.P_k = P_k