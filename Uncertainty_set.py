from gurobipy import *
from Assets import Asset_cls, create_assets
import random
import numpy as np
import csv
import math
import pandas as pd
random.seed(5)
np.random.seed(5)
# this class assumes uncertain parameters are equal through assets and time periods
class Uncertainty_set:
    #uncertain_parameters =pd.read_csv('Data/uncertain_parameters.csv').to_dict('records')[0]
    uncertain_parameters = {}
    with open('Data/uncertain_parameters.csv', mode='r') as inp:
        reader = csv.reader(inp)
        uncertain_parameters = {rows[0]: float(rows[1]) for rows in reader}
    def __init__(self,uncertain_parameters):
        self.d_bar = uncertain_parameters['Inherent_mean']
        self.d_hat = uncertain_parameters['Inherent_variation']
        self.zeta_bar = uncertain_parameters['Load_contribution_mean']
        self.zeta_hat = uncertain_parameters['Load_contribution_variation']
        self.gamma_bar = uncertain_parameters['C2C_contribution_mean']
        self.gamma_hat = uncertain_parameters['C2C_contribution_variation']

    @property
    def all_attributes(self):  # print all objects attributes as well as defined properties values
        return dict(vars(self))

    def get_gamma_sum(self,assets):
        #This function returns d(i,t), zeta(i,t) and gamma(j,i) for both bar and varaition using parameters above.
        gamma_sum = {}
        for i in range(len(assets)):
            if assets[i].pair!=-1:
                gamma_sum[i] = self.gamma_bar/self.gamma_hat
            else:
                gamma_sum[i] = 0
        return gamma_sum
    def get_parameters(self,number_of_assets,number_of_periods):
        #This function returns d(i,t), zeta(i,t) and gamma(j,i) for both bar and varaition using parameters above.
        #It assumes parameters are constant through time and machines
        parameters = {}
        parameters['d_bar'] = {(i, t): self.d_bar for i in range(number_of_assets) for t in
                                         range(number_of_periods)}
        parameters['d_hat'] = {(i, t): self.d_hat for i in range(number_of_assets) for t in
                                         range(number_of_periods)}
        parameters['zeta_bar'] = {(i, t): self.zeta_bar for i in range(number_of_assets) for t in
                                         range(number_of_periods)}
        parameters['zeta_hat'] = {(i, t): self.zeta_hat for i in range(number_of_assets) for t in
                                         range(number_of_periods)}
        parameters['gamma_bar'] = {(j, i): self.gamma_bar for j in range(number_of_assets) for i in
                                         range(number_of_assets) if j!=i}
        parameters['gamma_hat'] = {(j, i): self.gamma_hat for j in range(number_of_assets) for i in
                                         range(number_of_assets) if j!=i}
        return parameters
    def get_delta_prime(self,assets,budget,number_of_periods):
        delta = {(i, t): math.sqrt(t)*budget for i in range(len(assets)) for t in range(number_of_periods)}
        delta_prime = {}
        gamma_sum = self.get_gamma_sum(assets)
        u = self.get_parameters(len(assets),number_of_periods)
        for i in range(len(assets)):
            for t in range(number_of_periods):
                delta_prime[i,t] = delta[i, t]
                for tau in range(t):
                    delta_prime[i,t] += (u['d_bar'][i, tau] / u['d_hat'][i, tau])
                    delta_prime[i,t] += (u['zeta_bar'][i,tau] / u['zeta_hat'][i,tau])
                    delta_prime[i,t] += gamma_sum[i]
        return delta_prime
def create_unceratinty_set():
    return Uncertainty_set(Uncertainty_set.uncertain_parameters)
