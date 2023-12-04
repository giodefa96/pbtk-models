from Animals.Animal import Animal


class Trout(Animal):
    def __init__(self, B_w=1.0, Q_g=7.2, Q_c=2.07, V_cf= 0.098,
                 V_lc=0.012, V_sc=0.818, V_rc=0.063, V_kc=0.009,
                 Q_fc=0.085, Q_lc=0.029, Q_sc=0.600, Q_rc=0.230, Q_kc=0.056):
        
        # Parameters
        # Flows
        self.Q_f = Q_fc*Q_c
        self.Q_s = Q_sc*Q_c
        self.Q_r = Q_rc*Q_c
        self.Q_l = Q_lc*Q_c
        self.Q_k = Q_kc*Q_c
        
        # Volumes
        self.V_r = V_rc*B_w
        self.V_s = V_sc*B_w
        self.V_f = V_cf*B_w
        self.V_l = V_lc*B_w
        self.V_k = V_kc*B_w
        
        # Parameters
        self.Q_c = Q_c
        self.Q_g = Q_g
        
        # Tissue partition concentrations  
        self.C_a = 0
        self.C_k = 0
        self.C_l = 0  
        self.C_f = 0  
        self.C_s = 0  
        self.C_r = 0  
        self.C_v = 0
        self.C_m = 0 
