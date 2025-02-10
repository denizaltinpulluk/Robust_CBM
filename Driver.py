from Optimization import *
from Models import *
from Uncertainty_set import *
from Scenarios import *
from Assets import Asset_cls, create_assets

import time

start_time = time.time()

# Set random seed for reproducibility
random.seed(5)
np.random.seed(5)

# Initialize optimization and assets
opt = create_optimization_module()
assets = create_assets(opt)
U = create_uncertainty_set()

# Define budget levels
budget = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

# Loop through different budget levels
for i, b in enumerate(budget):
    scenarios = create_scenarios(number_of_scenarios=100, optimization=opt,
                                 uncertainty_set=U, assets=assets)
    # Set and solve robust optimization model
    opt.set_robust_model(assets=assets, budget = budget[i], Uncertainty_set=U)
    opt.solve_robust_model()

    # Write results to file
    results = open("results/results_%s_%.2f.txt" % ('comprehensive',budget[i]), "a")
    results.write("robust obj :\t %f\n" % (opt.robust_model.objVal))
    results.write("scenarios :\n")


    for scenario in scenarios:

        obj,  number_of_failure,total_penalty, = opt.solve_sim(assets,scenario)

        results.write("%s :\t %f\t %f\t %d\n" % (scenario.scenario_name, obj, total_penalty, number_of_failure))

    results.close()

    print()

end_time = time.time()

elapsed_time = end_time - start_time
print(f"Execution Time: {elapsed_time} seconds")
print()