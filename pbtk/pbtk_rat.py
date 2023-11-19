import matplotlib.pyplot as plt  
  
from scipy.integrate import odeint

from pbtk import AbstractPBTKModel

class RatPBTKModel(AbstractPBTKModel):  
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

    def __init__(self, Q_c=5.25, Q_p=5.25, Q_f=0.47, Q_l=1.31, Q_s=0.79, Q_r=2.68,   
                 V_f=0.022, V_l=0.012, V_s=0.174, V_r=0.012,   
                 P_f=56.72, P_l=4.64, P_s=1.54, P_r=4.64, P_b=18,   
                 Vmax=1.66, Km=0.55, C_inh_n=0.75,   
                 C_l=0, C_f=0, C_s=0, C_r=0, C_v=0,
                 integration_step=0.005,
                 exposure_time=4,
                 simulation_time=7):    
        
        # Parameters    
        self.Q_c, self.Q_p, self.Q_f, self.Q_l, self.Q_s, self.Q_r = Q_c, Q_p, Q_f, Q_l, Q_s, Q_r    
        self.V_f, self.V_l, self.V_s, self.V_r = V_f, V_l, V_s, V_r    
        self.P_f, self.P_l, self.P_s, self.P_r, self.P_b = P_f, P_l, P_s, P_r, P_b    
        self.Vmax, self.Km = Vmax, Km    
        self.integration_step = integration_step
        self.exposure_time = exposure_time
        self.simulation_time = simulation_time
    
        # Time points    
        self.t = [0.01 + self.integration_step*i for i in range(int(self.simulation_time/self.integration_step))]    
    
        # Initial conditions    
        self.C_inh_n = C_inh_n    
        self.C_l, self.C_f, self.C_s, self.C_r, self.C_v = C_l, C_f, C_s, C_r, C_v    
    
        # Lists to store results    
        self.C_a_list, self.C_v_list, self.C_l_list, self.C_f_list, self.C_s_list, self.C_r_list = [], [], [], [], [], []    
  
    def calculate_concentrations(self):  
        for time in self.t:  
            if time > self.exposure_time:  
                self.C_inh_n = 0  
  
            # Arterial blood  
            C_a = round((self.Q_c*self.C_v + self.Q_p*self.C_inh_n) / (self.Q_c + (self.Q_p/self.P_b)),2)  
  
            # Liver  
            dc_l_dt = round( (self.Q_l/self.V_l) * (C_a - self.C_l/self.P_l) - ((self.Vmax/self.V_l * self.C_l/self.P_l)/(self.Km + self.C_l/self.P_l)),2)  
            self.C_l = round(dc_l_dt * self.integration_step + self.C_l,2)  
  
            # Fat  
            dc_f_dt = round(self.Q_f/self.V_f * (C_a - self.C_f / self.P_f),2)  
            self.C_f = round(dc_f_dt * self.integration_step + self.C_f,2)  
  
            # Poorly perfused tissue  
            dc_s_dt = round(self.Q_s/self.V_s * (C_a - self.C_s / self.P_s),2)  
            self.C_s = round(dc_s_dt * self.integration_step + self.C_s,2)  
  
            # Richly perfused tissue  
            dc_r_dt = round(self.Q_r/self.V_r * (C_a - self.C_r / self.P_r),2)  
            self.C_r = round(dc_r_dt * self.integration_step + self.C_r,2)  
  
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
    model = RatPBTKModel()  
    model.calculate_concentrations()  
    model.save_plots('plots/plot1.png', 'plots/plot2.png')  