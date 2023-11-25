from abc import ABC, abstractmethod  
  
class AbstractPBTKModel(ABC):  
    @abstractmethod  
    def __init__(self):  
        pass  
  
    @abstractmethod  
    def calculate_concentrations(self):  
        pass  
  
    @abstractmethod  
    def plot_results(self):  
        pass  
  
    @abstractmethod  
    def save_plots(self, filename1, filename2):  
        pass  
    
    def calculate_integral(self, derivate, t, C):
        return derivate * t + C 