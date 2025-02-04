from gurobipy import *
from Optimization import *
import random
import numpy as np
import sys
import csv
import time
import math

penalty_factor = 10
random.seed(5)
np.random.seed(5)


def check_failure(optimization, assets, production, mp_r,scenario):
    deg = {}
    threshold = 100  # indicates what time asset exceed failure threshold
    omega = {}
    availability = {}
    M_time_left = {}
    PM_start = {}
    CM_start = {}
    fail = np.zeros(optimization.n_of_assets,dtype = bool)
    for i in range(optimization.n_of_assets):
        deg[i, 0] = assets[i].initial_degradation
        omega[i, 0] = assets[i].initial_degradation
        M_time_left[i] = 0
        fail[i] = False
        availability[i, 0] = True
        PM_start[i, 0] = 0
        CM_start[i, 0] = 0

    for t in range(optimization.n_of_weeks - 1):
        for i in range(optimization.n_of_assets):
            if M_time_left[i] == 0:
                deg[i, t + 1] = deg[i,t] + scenario.d[i, t] + scenario.zeta[i, t] * (
                        production[i, t].x / assets[i].Production_max) \
                                + scenario.gamma[assets[i].pair, i] * deg[assets[i].pair, t]
                availability[i, t + 1] = True
                PM_start[i, t + 1] = 0
                CM_start[i, t + 1] = 0
                if round(mp_r[i, t + 1].x) == 1:
                    M_time_left[i] = 2
                    deg[i, t + 1] = 0
                    availability[i, t + 1] = False
                    PM_start[i,t+1] = 1
                    CM_start[i, t + 1] = 0
                if deg[i, t + 1] > threshold:
                    fail[i] = True
                    M_time_left[i] = optimization.CM_Duration - 1
                    deg[i, t + 1] = 0
                    availability[i, t + 1] = False
                    CM_start[i, t + 1] = 1
                    PM_start[i, t + 1] = 0
            else:
                deg[i, t + 1] = 0
                M_time_left[i] = M_time_left[i] - 1
                availability[i, t + 1] = False
                PM_start[i, t + 1] = 0
                CM_start[i, t + 1] = 0

    return fail, availability, PM_start,CM_start

class Models:
    def create_robust_model(self, optimization, assets, budget, uncertainty_set):
        m_p, L, x, psi, m_c, u_f ,u_m, production, omega, omega_prime, p_prime = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
        v, z, v_zero, L_init, L_prime = {}, {}, {}, {}, {}
        pi_1, pi_2, pi_3, pi_4, pi_5, pi_6, pi_7 = {}, {}, {}, {}, {}, {}, {}
        model = Model("Robust_Model")
        u = uncertainty_set.get_parameters(optimization.n_of_assets, optimization.n_of_weeks)
        delta_prime = uncertainty_set.get_delta_prime(assets=assets, budget=budget,
                                                      number_of_periods=optimization.n_of_weeks)
        demand = optimization.demand
        print()
        for i in range(optimization.n_of_assets):
            L_init[i, 1] = assets[i].initial_degradation
            for k in range(2, optimization.n_of_maintenance):
                L_init[i, k] = 0

        for i in range(optimization.n_of_assets):
            for t in range(optimization.n_of_weeks):
                m_p[i, t] = model.addVar(
                    obj=optimization.preventive_maintenance * (1 - t * (1 / optimization.avg_lifetime)),
                    name="m_p_%d_%d" % (i, t), vtype=GRB.BINARY)
                # ------------------------------------------------------------------------------------------------------
                v_zero[i, t] = model.addVar(obj=0, name="v_zero_%d_%d" % (i, t), vtype=GRB.BINARY)

                # ------------------------------------------------------------------------------------------------------
                m_c[i, t] = model.addVar(obj=optimization.corrective_maintenance, name="m_c_%d_%d" % (i, t),
                                       vtype=GRB.BINARY)
                # m_c[i,t].start = maintenance['f'][i,t]
                # ------------------------------------------------------------------------------------------------------
                production[i, t] = model.addVar(obj=assets[i].production_cost, name="y_%d_%d" % (i, t),
                                                vtype=GRB.CONTINUOUS)
                # ------------------------------------------------------------------------------------------------------
                u_m[i, t] = model.addVar(obj=0, name="u_m_%d_%d" % (i, t), vtype=GRB.BINARY)
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for t in range(optimization.n_of_weeks):
                for k in range(1, optimization.n_of_maintenance):
                    v[i, t, k] = model.addVar(obj=0, name="v_%d_%d_%d" % (i, t, k), vtype=GRB.BINARY)
                    # v[i, t].start = maintenance['v'][i, t]
                    # ------------------------------------------------------------------------------------------------------
                    z[i, t, k] = model.addVar(obj=0, name="z_%d_%d_%d" % (i, t, k), vtype=GRB.BINARY)
                    # z[i, t].start = maintenance['z'][i, t]
                    # ------------------------------------------------------------------------------------------------------
                    L[i, t, k] = model.addVar(obj=0, name="L_Prime_%d_%d_%d" % (i, t, k), vtype=GRB.CONTINUOUS)
                    # ------------------------------------------------------------------------------------------------------
                    omega_prime[i, t, k] = model.addVar(obj=0, name="omega_prime_%d_%d_%d" % (i, t, k),
                                                        vtype=GRB.CONTINUOUS)
                    # ------------------------------------------------------------------------------------------------------
                    p_prime[i, t, k] = model.addVar(obj=0, name="p_prime_%d_%d_%d" % (i, t, k), vtype=GRB.CONTINUOUS)
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for t in range(optimization.n_of_weeks - 1):
                omega[i, t] = model.addVar(obj=0, name="omega_%d_%d" % (i, t), vtype=GRB.CONTINUOUS)
                u_f[i, t] = model.addVar(obj=0, name="u_f_%d_%d" % (i, t), vtype=GRB.BINARY)
                # u_f[i, t].start = maintenance['u_f'][i, t]
            u_f[i, optimization.n_of_weeks - 1] = model.addVar(obj=optimization.corrective_maintenance,
                                                                   name="u_f_%d_%d" % (i, optimization.n_of_weeks - 1),
                                                                   vtype=GRB.BINARY)
            # u_f[i, T - 1].start = maintenance['u_f'][i, T - 1]
            omega[i, optimization.n_of_weeks - 1] = model.addVar(obj=0.0001,
                                                                 name="omega_%d_%d" % (i, optimization.n_of_weeks - 1),
                                                                 vtype=GRB.CONTINUOUS)
        # ------------------------------------------------------------------------------------------------------
        for t in range(optimization.n_of_weeks):
            psi[t] = model.addVar(obj=assets[i].production_cost * penalty_factor, name="psi_%d" % (t), vtype=GRB.CONTINUOUS)
        # ------------------------------------------------------------------------------------------------------
        # Duals
        for i in range(optimization.n_of_assets):
            for t in range(1, optimization.n_of_weeks):
                for k in range(1, optimization.n_of_maintenance):
                    pi_1[i, t, k] = model.addVar(obj=0, name="pi_1_%d_%d_%d" % (i, t, k), vtype=GRB.CONTINUOUS)
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for t in range(1, optimization.n_of_weeks):
                for tau in range(1, t + 1):
                    for k in range(1, optimization.n_of_maintenance):
                        pi_2[i, t, tau, k] = model.addVar(obj=0, name="pi_2_%d_%d_%d_%d" % (i, t, tau, k),
                                                          vtype=GRB.CONTINUOUS)
                        pi_3[i, t, tau, k] = model.addVar(obj=0, name="pi_3_%d_%d_%d_%d" % (i, t, tau, k),
                                                          vtype=GRB.CONTINUOUS)
                        pi_4[i, t, tau, k] = model.addVar(obj=0, name="pi_4_%d_%d_%d_%d" % (i, t, tau, k),
                                                              vtype=GRB.CONTINUOUS)
                        pi_5[i, t, tau, k] = model.addVar(obj=0, name="pi_5_%d_%d_%d_%d" % (i, t, tau, k),
                                                              vtype=GRB.CONTINUOUS)
                        if assets[i].pair != -1:
                            pi_7[assets[i].pair, i, t, tau, k] = model.addVar(obj=0, name="pi_7_%d_%d_%d_%d_%d" %
                                                                                          (assets[i].pair, i, t, tau, k),
                                                                              vtype=GRB.CONTINUOUS)
                            pi_6[assets[i].pair, i, t, tau, k] = model.addVar(obj=0,
                                                                                 name="pi_6_%d_%d_%d_%d_%d" %
                                                                                      (assets[i].pair, i, t, tau, k),
                                                                                 vtype=GRB.CONTINUOUS)
        # ------------------------------------------------------------------------------------------------------
        # ------------------------------------------------------------------------------------------------------
        # Constraints
        for i in range(optimization.n_of_assets):
            for t in range(1, optimization.n_of_weeks):
                for tau in range(1, t + 1):
                    for k in range(1, optimization.n_of_maintenance):
                        model.addConstr((1 / u['d_hat'][i, t]) * pi_1[i, t, k] + pi_2[i, t, tau, k] - pi_3[i, t, tau, k]
                                        >= 1 - z[i, tau, k] - u_f[i, tau])
                        # ------------------------------------------------------------------------------------------------------
                        model.addConstr(
                                (1 / u['zeta_hat'][i, t]) * pi_1[i, t, k] + pi_4[i, t, tau, k] - pi_5[i, t, tau, k]
                                >= p_prime[i, tau, k] / assets[i].Production_max)
                        # ------------------------------------------------------------------------------------------------------
                        if assets[i].pair != -1:
                            model.addConstr(
                                (1 / u['gamma_hat'][assets[i].pair, i]) * pi_1[i, t, k] + pi_6[
                                    assets[i].pair, i, t, tau, k]
                                - pi_7[assets[i].pair, i, t, tau, k] >= omega_prime[assets[i].pair, tau, k])
                        # ------------------------------------------------------------------------------------------------------

        for i in range(optimization.n_of_assets):
            for t in range(1, optimization.n_of_weeks):
                for k in range(1, optimization.n_of_maintenance):
                    if assets[i].pair != -1:
                        model.addConstr(L[i, t, k] - L_init[i, k] >= delta_prime[i, t] * pi_1[i, t, k]
                                        + quicksum((u['d_bar'][i, t] + u['d_hat'][i, t]) * pi_2[i, t, tau, k]
                                                   - (u['d_bar'][i, t] )* pi_3[i, t, tau, k]
                                                   +(u['zeta_bar'][i, t] + u['zeta_hat'][i, t]) * pi_4[i, t, tau, k]
                                                   - (u['zeta_bar'][i, t])* pi_5[i, t, tau, k]
                                                   + (u['gamma_bar'][assets[i].pair, i] + u['gamma_hat'][assets[i].pair, i]) * pi_6[assets[i].pair, i, t, tau, k]
                                                   - (u['gamma_bar'][assets[i].pair, i] ) * pi_7[assets[i].pair, i, t, tau, k]
                                                   for tau in range(1, t + 1)),
                                        name='deg_prime_%d_%d_%d' % (i, t, k))
                    else:
                        model.addConstr(L[i, t, k] - L_init[i, k] >= delta_prime[i, t] * pi_1[i, t, k]
                                        + quicksum((u['d_bar'][i, t] + u['d_hat'][i, t]) * pi_2[i, t, tau, k]
                                                   - (u['d_bar'][i, t]) * pi_3[i, t, tau, k] +
                                                   (u['zeta_bar'][i, t] + u['zeta_hat'][i, t]) * pi_4[i, t, tau, k]
                                                   - (u['zeta_bar'][i, t]) * pi_5[i, t, tau, k]
                                                   for tau in range(1, t + 1)),
                                        name='deg_prime_%d_%d_%d' % (i, t, k))

        # End of Duality
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            model.addConstr(L[i, 0, 1] - L_init[i, 1] >= 0)
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for t in range(optimization.n_of_weeks):
                model.addConstr(quicksum(v[i, t, k] for k in range(1, optimization.n_of_maintenance)) <= 1)
                # ------------------------------------------------------------------------------------------------------
                model.addConstr(100 * u_f[i, t] <= omega[i, t], name='failure_con_%d_%d' % (i, t))
                # ------------------------------------------------------------------------------------------------------
                model.addConstr(1 - u_f[i, t] >= m_p[i, t])
                # ------------------------------------------------------------------------------------------------------
                model.addConstr(production[i, t] <= assets[i].Production_max * (1 - u_m[i, t]))
                # ------------------------------------------------------------------------------------------------------
                model.addConstr(production[i, t] <= assets[i].Production_max * (1 - u_f[i, t]))
                # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for t in range(optimization.n_of_weeks):
                for k in range(1, optimization.n_of_maintenance):
                    model.addConstr(v[i, t, k] <= m_p[i, t]+ m_c[i,t])
                    model.addConstr(
                        p_prime[i, t, k] >= production[i, t] - optimization.BigM2 * z[i, t, k] - optimization.BigM2 *
                        u_f[i, t],
                        name='p_prime_con1_%d_%d_%d' % (i, t, k))
                    model.addConstr(
                        omega[i, t] >= L[i, t, k] - optimization.BigM * z[i, t, k],
                        name='w_con1_%d_%d_%d' % (i, t, k))
                    if assets[i].pair != -1:
                        model.addConstr(
                            omega_prime[assets[i].pair, t, k] >= omega[assets[i].pair, t] - optimization.BigM * z[
                                i, t, k] - optimization.BigM * u_f[i, t],
                            name='w_prime_con1_%d_%d_%d' % (i, t, k))
                    model.addConstr(L[i, t, k] <= 100, name='threshold_%d_%d_%d' % (i, t, k))
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for k in range(1, optimization.n_of_maintenance):
                model.addConstr(quicksum(v[i, t, k] for t in range(optimization.n_of_weeks)) + v_zero[i, k] == 1)
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for k in range(2, optimization.n_of_maintenance):
                model.addConstr(v_zero[i, k] >= v_zero[i, k - 1])
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for t in range(optimization.n_of_weeks):
                model.addConstr(z[i, t, 1] <= quicksum(v[i, tau, 1] for tau in range(t + 1)) + u_m[i, t],
                                name='v_z_coupling_%d_%d_%d' % (i, t, 1))
                for k in range(2, optimization.n_of_maintenance):
                    model.addConstr(z[i, t, k] <= 1 - quicksum(v[i, tau, k - 1] for tau in range(t + 1))
                                    + quicksum(v[i, tau, k] for tau in range(t + 1)) + u_m[i, t],
                                    name='v_z_coupling_%d_%d_%d' % (i, t, k))
        # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            for t in range(1, optimization.n_of_weeks):
                model.addConstr(u_f[i, t - 1] - u_f[i, t] <= m_c[i, t])
        for t in range(optimization.n_of_weeks):
            model.addConstr(quicksum(u_m[i, t] for i in range(optimization.n_of_assets)) <= optimization.M_max)
            model.addConstr(
                quicksum(production[i, t] for i in range(optimization.n_of_assets)) >= demand[t] - psi[t],
                name="Demand_%d" % (t))
            # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            model.addConstr(m_p[i, 0] == 0)
            model.addConstr(m_c[i, 0] == 0)
            for t in range(optimization.n_of_weeks - optimization.CM_Duration, optimization.n_of_weeks):
                model.addConstr(m_p[i, t] == 0)
                model.addConstr(m_c[i, t] == 0)
            # ------------------------------------------------------------------------------------------------------
        for i in range(optimization.n_of_assets):
            model.addConstr(u_m[i, 0] == m_p[i, 0] + m_c[i, 0])
            model.addConstr(u_m[i, 1] == m_p[i, 1] + m_c[i, 1])
            model.addConstr(u_m[i, 2] == m_p[i, 2] + m_p[i, 1] + m_c[i, 2] + m_c[i, 1])
            model.addConstr(u_m[i, 3] == m_p[i, 3] + m_p[i, 2] + m_p[i, 1] + m_c[i, 3] + m_c[i, 2] + m_c[i, 1])
            model.addConstr(u_m[i, 4] == m_p[i, 4] + m_p[i, 3] + m_p[i, 2] + m_c[i, 4] + m_c[i, 3] + m_c[i, 2] + m_c[i, 1])
            for t in range(5, optimization.n_of_weeks):
                model.addConstr(
                    u_m[i, t] == m_p[i, t] + m_p[i, t - 1] + m_p[i, t - 2] + m_c[i, t] + m_c[i, t - 1] + m_c[i, t - 2] + m_c[
                        i, t - 3] + m_c[i, t - 4])
        #model.optimize()
        return model

    def solve_simulation(self, optimization, assets, scenario):
        robust_model = optimization.robust_model
        robust_model.update()
        asset_init_deg = []
        demand = optimization.demand
        PM_cost  = np.zeros(optimization.n_of_weeks)
        CM_cost = np.zeros(optimization.n_of_weeks)
        Prod_cost = np.zeros(optimization.n_of_weeks)
        penalty = np.zeros(optimization.n_of_weeks)
        total_cost = np.zeros(optimization.n_of_weeks)
        for i in range(optimization.n_of_assets):
            asset_init_deg.append( assets[i].initial_degradation)
        mp_r, production_r = {}, {}
        for i in range(optimization.n_of_assets):
            for t in range(optimization.n_of_weeks):
                mp_r[i, t] = robust_model.getVarByName("m_p_%d_%d" % (i, t))
                production_r[i, t] = robust_model.getVarByName("y_%d_%d" % (i, t))



        fail , availability,PM_start,CM_start = check_failure(optimization,assets, production_r,mp_r,scenario)


        for t in range(optimization.n_of_weeks):
            capacity = 0
            for i in range(optimization.n_of_assets):
                capacity += assets[i].Production_max * availability[i,t]
                PM_cost[t] += PM_start[i,t] * optimization.preventive_maintenance * (1 - t * (1 / optimization.avg_lifetime))
                CM_cost[t] += CM_start[i, t] * optimization.corrective_maintenance
                Prod_cost[t] += production_r[i, t].x * assets[i].production_cost
            if demand[t]<= capacity:
                penalty[t] = 0
            else:
                penalty[t] = (demand[t] - capacity)* assets[i].production_cost * penalty_factor
            total_cost[t] = PM_cost[t] + CM_cost[t] + Prod_cost[t] + penalty[t]
        obj = np.sum(total_cost)
        num_failure = sum(fail)
        total_penalty = np.sum(penalty)

        return obj,num_failure,total_penalty
