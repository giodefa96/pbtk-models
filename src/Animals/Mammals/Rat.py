
from Animals.Animal import Animal  
  
class Rat(Animal):  
    """    
    This class represents a rat in a Physiologically Based Toxicokinetic (PBTK) model.    
        
    The class inherits from the Animal class and initializes several physiological parameters     
    such as cardiac output (Q_c), alveolar ventilation rate (Q_p), blood flow to various tissues     
    (Q_f, Q_l, Q_s, Q_r), volume of various tissues (V_f, V_l, V_s, V_r), maximal velocity of     
    metabolism (V_max), and Michaelis–Menten affinity (K_m).    
        
    It also initializes tissue partition concentrations (C_l, C_f, C_s, C_r, C_v) to 0.    
        
    Parameters    
    ----------    
    Q_c : float, optional    
        Cardiac output, by default 5.25    
    Q_p : float, optional    
        Alveolar ventilation rate, by default 5.25    
    Q_f : float, optional    
        Fat blood flow, by default 0.47    
    Q_l : float, optional    
        Liver blood flow, by default 1.31    
    Q_s : float, optional    
        Poorly perfused tissues blood flow, by default 0.79    
    Q_r : float, optional    
        Richly perfused tissues blood flow, by default 2.68    
    V_f : float, optional    
        Fat volume, by default 0.022    
    V_l : float, optional    
        Liver volume, by default 0.012    
    V_s : float, optional    
        Poorly perfused tissues volume, by default 0.174    
    V_r : float, optional    
        Richly perfused tissues volume, by default 0.012    
    """    
    def __init__(self, Q_c=5.25, Q_p=5.25,  
                 Q_f=0.47, Q_l=1.31,  
                 Q_s=0.79, Q_r=2.68,  
                 V_f=0.022, V_l=0.012,  
                 V_s=0.174, V_r=0.012):  
        # Parameters  
        # Flows  
        self.Q_c = Q_c  
        self.Q_p = Q_p  
        self.Q_f = Q_f  
        self.Q_l = Q_l  
        self.Q_s = Q_s  
        self.Q_r = Q_r  
          
        # Volumes  
        self.V_f = V_f  
        self.V_l = V_l  
        self.V_s = V_s  
        self.V_r = V_r  
                  
        # Tissue partition concentrations  
        self.C_l = 0  
        self.C_f = 0  
        self.C_s = 0  
        self.C_r = 0  
        self.C_v = 0  
