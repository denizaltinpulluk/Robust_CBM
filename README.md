# Robust Condition-Based Fleet Operations & Maintenance

## Overview
This repository contains the implementation of a robust optimization model for operations and maintenance (O&M) decision-making. The model integrates multi-asset degradation rate interactions and operations-induced degradation to optimize fleet-wide O&M strategies. The optimized decisions are then evaluated under 100 different scenarios to assess performance. The outputs include O&M costs, total penalties, and the number of failures for each scenario.

## Installation & Requirements
- Python 3.x
- Gurobi (a valid Gurobi license is required to run the optimization model)
- Required Python libraries: (e.g., NumPy, Pandas, Matplotlib, etc.)

## Usage
### Running the Experiments
1. Ensure that the **results** folder exists and is empty before running the experiments.
2. Execute `Driver.py` to run all the experiments and generate results.
3. The code will produce output files in the `results` folder, with each file corresponding to different budget levels.

### Plotting Results
After obtaining the results:
1. Run `Experiments.py` to generate plots and analyze the outputs.

## File Structure
```
├── Data/                  # Contains input parameter files
│   ├── Asset_specs.csv
│   ├── model_parameters.csv
│   ├── uncertain_parameters.csv
│
├── results/               # Output directory (must be empty before running the code)
│   ├── <output_file_1>.txt
│
├── Driver.py              # Main script for running experiments
├── Experiments.py         # Script for analyzing and plotting results
├── Assets.py              # Defines asset properties and behavior
├── Optimization.py        # Implements optimization routines
├── Scenarios.py           # Generates and manages scenarios
├── Uncertainty_set.py     # Handles uncertainty modeling
├── Models.py              # Defines mathematical models for optimization
└── README.md              # Project documentation
```

## License
This code requires a valid **Gurobi license** to run the optimization model.

## Citation
If you use this code in your research, please consider citing the related paper (citation details will be provided once available).

---
For any questions or issues, feel free to open an issue on GitHub.

