
# H2O-H2 Rate Coefficients Calculator

This Python script calculates effective and thermal rate coefficients for rotational state-to-state transitions in H2O + H2 collisions using the method described by Daniel et al., with their data retrieved from BASECOL database. The same 27 values of temperature in the range between 5K and 1500K are considered as in the original work of Daniel et al. 

Note, in order to compute thermal rate coefficients for all transitions in ortho-H2O, the effective rate coefficients obtained for H2(j=2) projectile are used as an approximation for H2(j=4,6,8) projectiles. The same approach is used for some transitions in para-H2O, except those where the data for H2(j=2) projectile are missing. For those cases, the effective rate coefficients obtained for ortho-H2(j=1) projectile are used as an approximation for para-H2(j=2,4,6,8). This last feature was not available in the original code of Dubernet and coworkers.

## Acknowledgement

We acknowledge Dr. Marie Lise Dubernet for her significant contributions to this work.

## Citation

When using this code, please cite the following papers:

1. Daniel, Fabien, M-L. Dubernet, and Alain Grosjean. "Rotational excitation of 45 levels of ortho/para-H2O by excited ortho/para-H2 from 5 K to 1500 K: state-to-state, effective, and thermalized rate coefficients." Astronomy & Astrophysics 536 (2011): A76.

2. C. Joy, D. Bostan, B. Mandal and D. Babikov "Rate coefficients for rotational state-to-state transitions in H2O + H2 collisions as predicted by mixed quantum/classical theory" submitted to Astronomy & Astrophysics, 2024.

## Features

- Reads state-to-state rate coefficients from input files
- Computes thermal & effective rate coefficients for H2O-H2 collisions
- Handles both para and ortho forms of H2 and H2O
- Writes results to output files

## Requirements

- Python 3.x
- pandas
- numpy

## Usage

1. Ensure all required input files are in the `Dubernet_Data` directory.
2. Run the script:

   ```
   python compute_thermal_rates.py 
   ```

3. The script will generate output files in the `output-effective-rates` and `output-thermal-rates` directories.

## Input Files

The script requires the following input files in the `Dubernet_Data` directory:

- State-to-state rate coefficient files:
  - `H2O-para_H2-para_statetostate_ID148_v1.kij`
  - `H2O-para_H2-ortho_statetostate_ID149_v1.kij`
  - `H2O-ortho_H2-para_statetostate_ID94_v1.kij`
  - `H2O-ortho_H2-ortho_statetostate_ID147_v1.kij`

- Energy level files:
  - `H2-para_IDT4_vt1_IDC148_vc1.lev`
  - `H2-ortho_IDT16_vt1_IDC149_vc1.lev`
  - `H2O-para_IDT14_vt1_IDC149_vc1.lev`
  - `H2O-ortho_IDT1_vt1_IDC94_vc1.lev`

## Output Files

The script generates the following output files:

- Effective rates:
  - `effective_rates_para_h2o_para_h2.dat`
  - `effective_rates_para_h2o_ortho_h2.dat`
  - `effective_rates_ortho_h2o_para_h2.dat`
  - `effective_rates_ortho_h2o_ortho_h2.dat`

- Thermal rates:
  - `thermal_rates_output_para_h2o_para_h2.dat`
  - `thermal_rates_output_para_h2o_ortho_h2.dat`
  - `thermal_rates_output_ortho_h2o_para_h2.dat`
  - `thermal_rates_output_ortho_h2o_ortho_h2.dat`

## Functions

- `read_state_to_state_data`: Reads state-to-state rate coefficients from input files
- `read_energy_levels_hydrogen`: Reads energy levels of hydrogen from input files
- `read_energy_levels_water`: Reads energy levels of water from input files
- `compute_effective_rate_coefficients`: Computes effective rate coefficients
- `compute_thermal_rates`: Computes thermal rate coefficients with fallback logic
- `compute_thermal_rates_regular`: Computes thermal rate coefficients without cross-symmetry fallback
- `write_thermal_rates_to_file`: Writes thermal rate coefficients to output files
- `write_effective_rates_to_file`: Writes effective rate coefficients to output files

## Note

This script uses a fallback logic for missing rates in the para-H2 thermal rate calculation. For ortho-H2 and cross-symmetry calculations, it uses a regular method without fallback.
