import matplotlib.pyplot as plt  

from pbtk import AbstractPBTKModel
from Animals.Mammals.Rat import Rat

from Molecules.perchloroethylene import Perchloroethylene  

#export PYTHONPATH="${PYTHONPATH}:
class RatPBTKModel(AbstractPBTKModel, Rat):  
    """  
This class implements a four-compartment Physiologically Based Toxicokinetic (PBTK) model for a rat.   
  
The model includes the following compartments:   
- Arterial blood  
- Venous blood  
- Liver  
- Fat  
- Poorly perfused tissues (PPT)  
- Richly perfused tissues (RPT)  
  
The model uses the following parameters:  
- Cardiac output (Q_c)  
- Alveolar ventilation rate (Q_p)  
- Fat blood flow (Q_f)  
- Hepatic blood flow (Q_l)  
- Poorly perfused tissues (PPT) blood flow (Q_s)  
- Richly perfused tissues (RPT) blood flow (Q_r)  
- Fat volume (V_f)  
- Liver volume (V_l)  
- PPT volume (V_s)  
- RPT volume (V_r)  
- Fat:blood partition coefficient (P_f)  
- Liver:blood partition coefficient (P_l)  
- Poorly perfused tissue:blood partition coefficient (P_s)  
- Richly perfused tissue:blood partition coefficient (P_r)  
- Blood:air partition coefficient (P_b)  
- Maximal velocity of metabolism (Vmax)  
- Michaelis–Menten affinity (Km)  
  
The model calculates the concentration of a substance in each compartment over time, given an initial concentration in the arterial blood. The concentration in each compartment is updated at each time step using the differential equations that describe the flow of the substance between compartments.  
  
The results of the model can be plotted using the `plot_results` method, which creates two plots: one showing the concentration over time in the arterial and venous blood, liver, PPT, and RPT, and another showing the concentration over time in the fat.  
  
The model can also save the plots to files using the `save_plots` method.  
"""  

    def __init__(self,
                 rat,
                 molecule,
                 C_inh_n=0.75,
                 integration_step=0.005,
                 exposure_time=4,
                 simulation_time=7): 
        
        
        super().__init__()  # Call the parent class constructor if necessary  
        self.rat = rat  # Initialize the Rat class  
  
        # Now you can access the parameters of the Rat class using the rat object  
        self.Q_c = self.rat.Q_c  
        self.Q_p = self.rat.Q_p  
        self.Q_f = self.rat.Q_f  
        self.Q_l = self.rat.Q_l  
        self.Q_s = self.rat.Q_s  
        self.Q_r = self.rat.Q_r  
        self.V_f = self.rat.V_f  
        self.V_l = self.rat.V_l  
        self.V_s = self.rat.V_s  
        self.V_r = self.rat.V_r  
        
        
        
        self.molecule = molecule  # Initialize the Molecule class  
        
        # Metabolism parameters  
        self.V_max, self.K_m = self.molecule.V_max, self.molecule.K_m
        
        # Parameters molecule interaction with tissues
        self.P_f, self.P_l, self.P_s = self.molecule.P_f, self.molecule.P_l, self.molecule.P_s
        self.P_r, self.P_b =  self.molecule.P_r, self.molecule.P_b   
         
        self.integration_step = integration_step
        self.exposure_time = exposure_time
        self.simulation_time = simulation_time
    
        # Time points    
        self.t = [0.01 + self.integration_step*i for i in range(int(self.simulation_time/self.integration_step))]    
    
        # Initial conditions    
        self.C_inh_n = C_inh_n    
        self.C_l, self.C_f, self.C_s, self.C_r, self.C_v = self.rat.C_l, self.rat.C_f, self.rat.C_s, self.rat.C_r, self.rat.C_v    
    
        # Lists to store results    
        self.C_a_list, self.C_v_list, self.C_l_list, self.C_f_list, self.C_s_list, self.C_r_list = [], [], [], [], [], []    
  
    def calculate_concentrations(self):  
        for time in self.t:  
            if time > self.exposure_time:  
                self.C_inh_n = 0  
  
            # Arterial blood  
            C_a = round((self.Q_c*self.C_v + self.Q_p*self.C_inh_n) / (self.Q_c + (self.Q_p/self.P_b)),2)  
  
            # Liver  
            dc_l_dt = round( (self.Q_l/self.V_l) * (C_a - self.C_l/self.P_l) - ((self.V_max/self.V_l * self.C_l/self.P_l)/(self.K_m + self.C_l/self.P_l)),2)  
            self.C_l = self.calculate_integral(dc_l_dt, self.integration_step, self.C_l)
  
            # Fat  
            dc_f_dt = round(self.Q_f/self.V_f * (C_a - self.C_f / self.P_f),2)  
            self.C_f = self.calculate_integral(dc_f_dt, self.integration_step, self.C_f)
  
            # Poorly perfused tissue  
            dc_s_dt = round(self.Q_s/self.V_s * (C_a - self.C_s / self.P_s),2)  
            self.C_s = self.calculate_integral(dc_s_dt, self.integration_step, self.C_s)
  
            # Richly perfused tissue  
            dc_r_dt = round(self.Q_r/self.V_r * (C_a - self.C_r / self.P_r),2)  
            self.C_r = self.calculate_integral(dc_r_dt, self.integration_step, self.C_r)
  
            # Venous blood  
            self.C_v = round((self.Q_l * self.C_l/self.P_l + self.Q_f*self.C_f/self.P_f + self.Q_s * self.C_s/self.P_s + self.Q_r * self.C_r/self.P_r) / (self.Q_c),2)  
  
            # Append results to lists  
            self.C_a_list.append(C_a),self.C_v_list.append(self.C_v),self.C_l_list.append(self.C_l)  
            self.C_f_list.append(self.C_f),self.C_s_list.append(self.C_s),self.C_r_list.append(self.C_r)  
  
        return self.C_a_list, self.C_v_list, self.C_l_list, self.C_f_list, self.C_s_list, self.C_r_list  
    
    def plot_results(self):  
        fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 12))  
        
        ax1.plot(self.t, self.C_v_list, label='Venous')  
        ax1.plot(self.t, self.C_a_list, label='Arterial blood')  
        ax1.plot(self.t, self.C_l_list, label='Liver')  
        ax1.plot(self.t, self.C_s_list, label='Poorly perfused tissue')
        ax1.plot(self.t, self.C_r_list, label='Richly perfused tissue')
        ax1.set_xlabel('Time (h)')  
        ax1.set_ylabel('Concentration (mg/L)')  
        ax1.set_title('Concentration over time in Venous, Arterial blood, Liver, Poorly perfused tissue and Richly perfused tissue')  
        ax1.legend()  
  
        ax2.plot(self.t, self.C_f_list, label='Fat')  
        ax2.set_xlabel('Time (h)')  
        ax2.set_ylabel('Concentration (mg/L)')  
        ax2.set_title('Concentration over time in Fat')  
        ax2.legend()  
  
        plt.tight_layout()  
        plt.show() 

    def save_plots(self, filename1, filename2):  
        # First plot  
        fig1, ax1 = plt.subplots(figsize=(10,6))  
        ax1.plot(self.t, self.C_v_list, label='Venous')    
        ax1.plot(self.t, self.C_a_list, label='Arterial blood')    
        ax1.plot(self.t, self.C_l_list, label='Liver')    
        ax1.plot(self.t, self.C_s_list, label='Poorly perfused tissue')  
        ax1.plot(self.t, self.C_r_list, label='Richly perfused tissue')  
        ax1.set_xlabel('Time (h)')    
        ax1.set_ylabel('Concentration (mg/L)')    
        ax1.set_title('Concentration over time in Venous, Arterial blood, Liver, Poorly perfused tissue and Richly perfused tissue')    
        ax1.legend()    
        fig1.savefig(filename1)  
        plt.close(fig1)  # Close the figure  
    
        # Second plot  
        fig2, ax2 = plt.subplots(figsize=(10,6))  
        ax2.plot(self.t, self.C_f_list, label='Fat')    
        ax2.set_xlabel('Time (h)')    
        ax2.set_ylabel('Concentration (mg/L)')    
        ax2.set_title('Concentration over time in Fat')    
        ax2.legend()    
        fig2.savefig(filename2)  
        plt.close(fig2)  # Close the figure  


if __name__ == '__main__':
    rat = Rat()
    perchloroethylene = Perchloroethylene()
    model = RatPBTKModel(rat, perchloroethylene)  
    model.calculate_concentrations()  
    model.plot_results()  