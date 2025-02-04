from Optimization import *
from Models import *
from Uncertainty_set import *
from Scenarios import *
from Assets import Asset_cls, create_assets
random.seed(5)
np.random.seed(5)



opt = create_optimization_module()
assets = create_assets(opt)
U = create_unceratinty_set()
budget = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
for i in range(9):
    scenarios = create_scenarios(number_of_scenarios=100, optimization=opt,
                                 uncertainty_set=U, assets=assets)
    opt.set_robust_model(assets=assets, budget = budget[i], Uncertainty_set=U)
    opt.solve_robust_model()

    results = open("results/results_%s_%.2f.txt" % ('comprehensive',budget[i]), "a")
    results.write("robust obj :\t %f\n" % (opt.robust_model.objVal))
    results.write("scenarios :\n")


    for scenario in scenarios:

        obj,  number_of_failure,total_penalty, = opt.solve_sim(assets,scenario)

        results.write("%s :\t %f\t %f\t %d\n" % (scenario.scenario_name, obj, total_penalty, number_of_failure))

    results.close()

    print()
print()