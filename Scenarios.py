from gurobipy import *

import random
import numpy as np

from Optimization import *
from Uncertainty_set import *

# Set random seed for reproducibility
random.seed(5)
np.random.seed(5)

def generate_scenario(optimization,uncertainty_set,assets,scenario_number):
    """Generates a single scenario based on the uncertainty set"""
    scenario = {'d': {}, 'zeta': {}, 'gamma': {},'number':scenario_number}
    for i in range(optimization.n_of_assets):
        for t in range(optimization.n_of_weeks):
            scenario['d'][i, t] = np.random.normal(uncertainty_set.d_bar, uncertainty_set.d_hat)
            scenario['zeta'][i, t] = np.random.normal(uncertainty_set.zeta_bar, uncertainty_set.zeta_hat)
        if assets[i].pair != -1:
            scenario['gamma'][assets[i].pair, i] = np.random.normal(uncertainty_set.gamma_bar,
                                                                       uncertainty_set.gamma_hat)
    return scenario
class Scenarios:
    def __init__(self,scenario):
        self.scenario_name = scenario['number']
        self.d = scenario['d']
        self.zeta = scenario['zeta']
        self.gamma = scenario['gamma']

    @property
    def all_attributes(self):  # print all objects attributes as well as defined properties values
        return dict(vars(self))
def create_scenarios(number_of_scenarios,optimization,uncertainty_set,assets):
    """Creates multiple scenarios based on input parameters."""
    s = []
    for i in range(number_of_scenarios):
        scenario = generate_scenario(optimization,uncertainty_set,assets,i+1)
        s.append(Scenarios(scenario))
    return s

