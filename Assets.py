import os
import glob
import pandas as pd
import numpy as np
import random

random.seed(5)

class Asset_cls:
    Assets_specifications=pd.read_csv("Data/Asset_specs.csv").T.to_dict()

    #***********************************************Objects attributes**************************************************
    def __init__(self,N,P_max,cost,pair):
        self.Name=N
        self.Production_max=P_max
        self.production_cost=cost
        self.initial_degradation =(5)+(75)*random.random()
        self.pair = pair

    @property
    def all_attributes(self):  # print all objects attributes as well as defined properties values
        return dict(vars(self))



    #**************************Setters and Getters**************************
def create_assets(opt):
    assets = []
    for v in range(0, opt.n_of_assets):
        if v % 2 == 0:
            pair = v + 1
        elif v % 2 == 1:
            pair= v - 1
        else:
            pair = -1 #no pair
        assets.append(Asset_cls("Asset_%d" % (v + 1),
                                       Asset_cls.Assets_specifications[v]['Max_Output'],
                                       Asset_cls.Assets_specifications[v]['Prod_Cost'],
                                       pair))  # create Asset objects(Name,age,profile)
    return assets
if __name__ == '__main__':
    Asset=create_assets()
    for i in Asset:
        print(i.all_attributes)