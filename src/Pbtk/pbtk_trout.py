import matplotlib.pyplot as plt  

from pbtk import AbstractPBTKModel
from Animals.Fishes.Trout import Trout
from Molecules.Tetrachloroethane import Tetrachloroethane


import matplotlib.pyplot as plt  


class TroutePBTKModel(AbstractPBTKModel):  
    """  
    A class used to represent a Physiologically Based Toxicokinetic (PBTK model for
    1,1,2,2-tetrachloroethane in rainbow trout.  
  
    ...  
  
    Attributes  
    ----------  
    Q_f, Q_s, etc. : float  
        Blood flow rates to different parts of the body.  
    V_r, V_s, etc. : float  
        Volumes of different tissues in the body.  
    Q_c, P_b, etc. : float  
        Other model parameters.  
    C_a_list, C_l_list, etc. : list  
        Lists to store the concentrations in different tissues over time.  
  
    Methods  
    -------  
    calculate_concentrations():  
        Calculates the concentrations in different tissues over time.  
    plot_results():  
        Creates a series of plots showing the concentrations in different tissues over time.  
    save_plots(filename):  
        Saves the plots to a file.  
    """  
    def __init__(self,
                 trout,
                 molecule,
                 T_stop=70, Dur_exp=48, C_int=0.25, C_ont=1.06,
                 integration_step=0.05):
        
        # Initialize the Trout class
        self.trout = trout
        
        # Parameters
        # Flows
        self.Q_f = self.trout.Q_f
        self.Q_s = self.trout.Q_s
        self.Q_r = self.trout.Q_r
        self.Q_l = self.trout.Q_l
        self.Q_k = self.trout.Q_k
        
        # Volumes
        self.V_r = self.trout.V_r
        self.V_s = self.trout.V_s
        self.V_f = self.trout.V_f
        self.V_l = self.trout.V_l
        self.V_k = self.trout.V_k
        
        # Parameters
        self.Q_c = self.trout.Q_c
        self.Q_g = self.trout.Q_g
        
        
        # Initialize the Molecule class
        self.molecule = molecule
        
        self.P_b = self.molecule.P_b
        self.P_f = self.molecule.P_f
        self.P_l = self.molecule.P_l
        self.P_r = self.molecule.P_r
        self.P_s = self.molecule.P_s
        self.P_k = self.molecule.P_k
        self.V_max = self.molecule.V_max
        self.K_m =self.molecule.K_m
        
        self.integration_step = integration_step
        self.exposure_time = Dur_exp
        self.simulation_time = T_stop
        
        # Time points    
        self.t = [0.05 + self.integration_step*i for i in range(int(self.simulation_time/self.integration_step))]
        # Initial conditions
        self.C_insp = C_ont
        self.C_int = C_int
        self.C_a, self.C_m = self.trout.C_a, self.trout.C_m
        self.C_f, self.C_s, self.C_r = self.trout.C_f, self.trout.C_s, self.trout.C_r
        self.C_l, self.C_k, self.C_v = self.trout.C_l, self.trout.C_k, self.trout.C_v
        
        # Lists to store results
        self.C_a_list, self.C_v_list, self.C_l_list, self.C_f_list, self.C_s_list, self.C_r_list, self.C_k_list = [], [], [], [], [], [], []
        
        
    def calculate_concentrations(self):
        for time in self.t:
            if time > self.exposure_time:
                self.C_insp = 0
            
            # Gill uptake limitation
            if self.Q_c*self.P_b > self.Q_r:
                GUL = self.Q_g
            else:
                GUL = self.Q_c*self.P_b
            
            # Model equationS
            # Calculation of blood concentration
            # Arterial blood concentration (mg/L)
            C_a = self.C_v+GUL*(self.C_insp-self.C_v/self.P_b)/self.Q_c
            # Venous blood concentration (mg/L)
            self.C_v = (self.Q_f*self.C_f/self.P_f + (self.Q_l + self.Q_r) * self.C_l/self.P_l + 0.4*self.Q_s*self.C_s/self.P_s + (self.Q_k + 0.6*self.Q_s) * self.C_k/self.P_k) / self.Q_c
            
            # Liver metabolism
            dc_met_dt = (self.V_max*self.C_l/self.P_l)/(self.K_m + self.C_l/self.P_l) / self.V_l
            self.C_m = self.calculate_integral(dc_met_dt, self.integration_step, self.C_m)
            # Liver
            dcl_dt = (self.Q_l*C_a + self.Q_r*self.C_r/self.P_r)/self.V_l - (self.Q_l+self.Q_r)/self.V_l * self.C_l/self.P_l - dc_met_dt
            self.C_l = self.calculate_integral(dcl_dt, self.integration_step, self.C_l)

            
            # Fat Tissue
            dc_f_dt = self.Q_f/self.V_f*(C_a-self.C_f/self.P_f)
            self.C_f = self.calculate_integral(dc_f_dt, self.integration_step, self.C_f)
            
            # Richly perfused tissue
            dc_r_dt = self.Q_r/self.V_r*(C_a-self.C_r/self.P_r)
            self.C_r = self.calculate_integral(dc_r_dt, self.integration_step, self.C_r)
            
            # Poorly perfused tissue
            dc_s_dt = self.Q_s/self.V_s*(C_a-self.C_s/self.P_s)
            self.C_s = self.calculate_integral(dc_s_dt, self.integration_step, self.C_s)
            
            # Kidney
            dc_k_dt = ((self.Q_k * C_a + 0.6 * self.Q_s * self.C_s/self.P_s) / self.V_k) - ((self.Q_k + 0.6 * self.Q_s) / self.V_k) * (self.C_k/self.P_k)
            self.C_k = self.calculate_integral(dc_k_dt, self.integration_step, self.C_k)
        
            # Append results to lists
            self.C_a_list.append(C_a)
            self.C_l_list.append(self.C_l)
            self.C_f_list.append(self.C_f)
            self.C_s_list.append(self.C_s)
            self.C_r_list.append(self.C_r)
            self.C_k_list.append(self.C_k)
            self.C_v_list.append(self.C_v)
            
        return self.C_a_list, self.C_l_list, self.C_f_list, self.C_s_list, self.C_r_list, self.C_k_list
    
        
    def plot_results(self):  
        fig, axs = plt.subplots(7, figsize=(10, 30))  
    
        axs[0].plot(self.t, self.C_a_list, label='C_a')  
        axs[0].set_title('Arterial blood concentration over time')  
        axs[0].set_xlabel('Time (h)')  
        axs[0].set_ylabel('Concentration (mg/L)')  
    
        axs[1].plot(self.t, self.C_l_list, label='C_l')  
        axs[1].set_title('Liver concentration over time')  
        axs[1].set_xlabel('Time (h)')  
        axs[1].set_ylabel('Concentration (mg/L)')  
    
        axs[2].plot(self.t, self.C_f_list, label='C_f')  
        axs[2].set_title('Fat Tissue concentration over time')  
        axs[2].set_xlabel('Time (h)')  
        axs[2].set_ylabel('Concentration (mg/L)')  
    
        axs[3].plot(self.t, self.C_s_list, label='C_s')  
        axs[3].set_title('Poorly perfused tissue concentration over time')  
        axs[3].set_xlabel('Time (h)')  
        axs[3].set_ylabel('Concentration (mg/L)')  
    
        axs[4].plot(self.t, self.C_r_list, label='C_r')  
        axs[4].set_title('Richly perfused tissue concentration over time')  
        axs[4].set_xlabel('Time (h)')  
        axs[4].set_ylabel('Concentration (mg/L)')  
    
        axs[5].plot(self.t, self.C_k_list, label='C_k')  
        axs[5].set_title('Kidney concentration over time')  
        axs[5].set_xlabel('Time (h)')  
        axs[5].set_ylabel('Concentration (mg/L)')  
    
        axs[6].plot(self.t, self.C_v_list, label='C_v')  
        axs[6].set_title('Venous blood concentration over time')  
        axs[6].set_xlabel('Time (h)')  
        axs[6].set_ylabel('Concentration (mg/L)')  
    
        for ax in axs.flat:  
            ax.legend()  
    
        plt.tight_layout()  
        plt.show()  

    
    def save_plots(self, filename):  
        fig, axs = plt.subplots(7, figsize=(10, 30))  
    
        axs[0].plot(self.t, self.C_a_list, label='C_a')  
        axs[0].set_title('Arterial blood concentration over time')  
        axs[0].set_xlabel('Time (h)')  
        axs[0].set_ylabel('Concentration (mg/L)')  
    
        axs[1].plot(self.t, self.C_l_list, label='C_l')  
        axs[1].set_title('Liver concentration over time')  
        axs[1].set_xlabel('Time (h)')  
        axs[1].set_ylabel('Concentration (mg/L)')  
    
        axs[2].plot(self.t, self.C_f_list, label='C_f')  
        axs[2].set_title('Fat Tissue concentration over time')  
        axs[2].set_xlabel('Time (h)')  
        axs[2].set_ylabel('Concentration (mg/L)')  
    
        axs[3].plot(self.t, self.C_s_list, label='C_s')  
        axs[3].set_title('Poorly perfused tissue concentration over time')  
        axs[3].set_xlabel('Time (h)')  
        axs[3].set_ylabel('Concentration (mg/L)')  
    
        axs[4].plot(self.t, self.C_r_list, label='C_r')  
        axs[4].set_title('Richly perfused tissue concentration over time')  
        axs[4].set_xlabel('Time (h)')  
        axs[4].set_ylabel('Concentration (mg/L)')  
    
        axs[5].plot(self.t, self.C_k_list, label='C_k')  
        axs[5].set_title('Kidney concentration over time')  
        axs[5].set_xlabel('Time (h)')  
        axs[5].set_ylabel('Concentration (mg/L)')  
    
        axs[6].plot(self.t, self.C_v_list, label='C_v')  
        axs[6].set_title('Venous blood concentration over time')  
        axs[6].set_xlabel('Time (h)')  
        axs[6].set_ylabel('Concentration (mg/L)')  
    
        for ax in axs.flat:  
            ax.legend()  
    
        plt.tight_layout()  
        plt.savefig(filename)  
        plt.close(fig)  
        

if __name__ == '__main__':
    # Instantiate the trout
    trout = Trout()
    # Instantiate the molecule
    molecule = Tetrachloroethane()
    # Instantiate the model
    model = TroutePBTKModel(trout=trout, molecule=molecule)
    
    # Run the model
    model.calculate_concentrations()
    
    # Plot the results
    model.plot_results()
    