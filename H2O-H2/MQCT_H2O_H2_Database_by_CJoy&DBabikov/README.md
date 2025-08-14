
## H2O-H2 Rate Coefficients Calculator

This Python script calculates effective and thermal rate coefficients for rotational state-to-state transitions in H2O + H2 collisions using mixed quantum-classical theory (MQCT) data. The calculations are performed for 21 temperature values ranging from 10K to 2000K.

Note, the effective rate coefficients are computed by summing over all final H2 rotational states (J2') for each transition. For thermal rate coefficients, a Boltzmann averaging procedure is applied over initial H2 rotational states.

## Citation

When using this code, please cite the following paper:

Joy, C., Bostan, D., Mandal, B., & Babikov, D. (2024). Rate coefficients for rotational state-to-state transitions in H2O+ H2 collisions as predicted by mixed quantum–classical theory. *Astronomy & Astrophysics*, 692, A229.

## Features

- Reads state-to-state rate coefficients from MQCT data files
- Computes effective rate coefficients by summing over final H2 states
- Computes thermal rate coefficients with Boltzmann averaging
- Handles ortho/para-H2O + ortho/para-H2 collision system
- Writes results to formatted output files

## Requirements

- Python 3.6+
- pandas
- numpy

## Usage

1. Ensure the input data file is in the `MQCT_DATA` directory.
2. Run the script:

   ```
   python Compute_Effective_thermal_rates.py
   ```

3. The script will generate output files in the `output-effective-rates` and `output-thermal-rates` directories.

## Input Files

The script requires the following input file in the `MQCT_DATA` directory:

- State-to-state rate coefficient file:
  - `MQCT_data_paraH2O_paraH2.txt`
  - Note: If one needs to compute for other symmetries like ortho-H2O + para-H2, ortho-H2O + ortho-H2, or para-H2O + ortho-H2, change the input file accordingly (e.g., MQCT_data_orthoH2O_paraH2.txt, etc.).

## Output Files

The script generates the following output files:

- Effective rates:
  - `output_effective_rates.dat` (in `output-effective-rates` directory)

- Thermal rates:
  - `output_thermal_rates.dat` (in `output-thermal-rates` directory)

## Functions

- `read_state_to_state_data`: Parses MQCT data file and extracts rate coefficients
- `compute_effective_rates`: Computes effective rate coefficients by summing over final H2 states
- `compute_thermal_rates`: Computes thermal rate coefficients with Boltzmann averaging and fallback logic
- `write_effective_rates`: Writes effective rate coefficients to formatted output file
- `write_thermal_rates`: Writes thermal rate coefficients to formatted output file

## Author

Carolin  
January 2025  
Version: 1.0
