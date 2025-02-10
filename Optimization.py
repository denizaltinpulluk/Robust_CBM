from gurobipy import *

from Models import *

import random
import numpy as np
import csv
from Uncertainty_set import *
from Scenarios import *
import time
import math
import pandas as pd

# Set random seed for reproducibility
random.seed(5)
np.random.seed(5)

def reader():
    parameters = {}
    with open('Data/model_parameters.csv', mode='r') as inp:
        reader = csv.reader(inp)
        parameters = {rows[0]: float(rows[1]) for rows in reader}
    return parameters

class Optimization:
    parameters = reader()
    def __init__(self,parameters):
        self.n_of_assets = int(parameters['Number_of_assets'])
        self.n_of_weeks = int(parameters['Number_of_weeks'])
        self.n_of_maintenance = int(parameters['Number_of_cycles'])
        self.MP_mipgap = parameters['Robust_MIPGAP']
        self.scenario_mipgap = parameters['Simulation_MIPGAP']
        self.time_limit = parameters['Simulation_time_limit']
        self.BigM = parameters['BigM_for_degradation']
        self.BigM2 = parameters['BigM_for_production']
        self.preventive_maintenance = parameters['Preventive_maintenance_cost']
        self.corrective_maintenance = parameters['Corrective_maintenance_cost']
        self.CM_Duration = int(parameters['Corrective_maintenance_duration'] ) # Corrective maintenance duration
        self.M_max = int(parameters['Maintenance_capacity'])
        self.avg_lifetime = int(parameters['Average_lifetime'])
        self.robust_model = None
        self.deterministic_model = None
        self.perfect_value_information_obj = []
        self.perfect_value_information_maintenance = []
        self.simulation_obj = []
        self.simulation_maintenance = []
        self.penalty =None
        self.failure =None
        self.demand = {(t): 800 + (150) * random.random() for t in range(int(parameters['Number_of_weeks']))}

    @property
    def all_attributes(self):  # print all objects attributes as well as defined properties values
        return dict(vars(self))

    def set_robust_model(self, assets, budget,Uncertainty_set):
        self.robust_model = Models().create_robust_model(self, assets, budget, Uncertainty_set)
        self.robust_model.update()

    def solve_robust_model(self):
        self.robust_model.setParam("MIPGap", self.MP_mipgap)
        self.robust_model.optimize()
        self.robust_model.update()

    def solve_sim(self,assets,scenario):
        obj,num_failure,total_penalty = Models().solve_simulation(self,assets,scenario)
        return obj, num_failure, total_penalty




def create_optimization_module():
    return Optimization(Optimization.parameters)


