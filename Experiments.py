import os
import numpy as np
import matplotlib.pyplot as plt


def plot_results(experiments):
    """Plots results from multiple experiment runs."""
    robus_objs = []
    sim_obj_means = []
    budget = []
    num_failure = []
    penalty = []
    for exp in experiments:
        budget.append(exp.budget)
        robus_objs.append(exp.robust_obj)
        sim_obj_means.append(np.mean(exp.simulation_objs))
        num_failure.append(np.mean(exp.simulation_failure))
        penalty.append(np.mean(exp.simulation_penalty))

    fig = plt.figure()
    plt.plot(budget, sim_obj_means, label=' Sim. Ave')
    plt.plot(budget, robus_objs, label=' robust')
    plt.xlabel("Budget")
    plt.ylabel("Obj F. Value")
    plt.legend()
    #plt.show()
    plt.savefig('OM_Costs.jpg')


    fig = plt.figure()
    plt.plot(budget, penalty, label=' Sim. Ave')
    plt.xlabel("Budget")
    plt.ylabel("Total Penalty")
    plt.legend()
    #plt.show()
    plt.savefig('penalty.jpg')


    fig = plt.figure()
    plt.plot(budget, num_failure, label=' Sim. Ave')
    plt.xlabel("Budget")
    plt.ylabel("Number of Failures")
    plt.legend()
    #plt.show()

    plt.savefig('num_failure.jpg')
class Experiments:
    """Stores results from experiment files and provides structured access to attributes."""
    def __init__(self,results):
        self.name = results['experiment_name']
        self.budget = results['budget']
        self.robust_obj = results['robust_obj']
        self.simulation_objs = results['simulations']
        self.simulation_penalty= results['penalty']
        self.simulation_failure = results['failure']



def create_experiments_from_files():
    """Reads results from text files in 'results' directory and creates Experiment objects."""
    path = "results"
    files = os.listdir(path)
    experiments = []

    for file in files:
        results = {}
        separator = file.split("_")
        results['experiment_name'] = separator[1]
        results['budget'] = float(separator[2][0:4])
        path_to_file = path + '/' + file
        simulation_objs = []
        simulation_failures =[]
        simulation_penalty = []
        with open(path_to_file) as f:
            contents = f.readlines()
        lines = []
        for content in contents:
            line = content.split()
            lines.append(line)
        results['robust_obj'] = float(lines[0][3])
        for i in range(2,len(lines)):
            simulation_objs.append(float(lines[i][2]))
            simulation_penalty.append(float(lines[i][3]))
            simulation_failures.append(int(lines[i][4]))
        results['simulations'] = np.array(simulation_objs)
        results['penalty'] = np.array(simulation_penalty)
        results['failure'] = np.array(simulation_failures)
        experiments.append(Experiments(results))
    return experiments
if __name__ == '__main__':
    experiments = create_experiments_from_files()
    plot_results(experiments)
    print()