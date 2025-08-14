"""
H2O-H2 Rate Coefficients Calculator

This script calculates effective and thermal rate coefficients for H2O-H2 collisions according to the prescription of Daniel et al using their data retrieved from BASECOL database.

Acknowledgement:
We acknowledge Dr. Marie Lise Dubernet for her significant contributions.

When using this code, please cite the following papers:

1. Daniel, Fabien, M-L. Dubernet, and Alain Grosjean. "Rotational excitation of 45 levels
   of ortho/para-H2O by excited ortho/para-H2 from 5 K to 1500 K: state-to-state, effective,
   and thermalized rate coefficients." Astronomy & Astrophysics 536 (2011): A76.

2. C. Joy, D. Bostan, B. Mandal and D. Babikov "Rate coefficients for rotational state-to-state
   transitions in H2O + H2 collisions as predicted by mixed quantum/classical theory"
   submitted to Astronomy & Astrophysics, 2024. 
"""

import pandas as pd 
import numpy as np

# Read the state-to-state rate coefficients from the file.
def read_state_to_state_data(filename):
    """
    Read state-to-state rate coefficients from a file.
    
    Args:
    filename (str): Path to the input file.
    
    Returns:
    tuple: Contains temperature array, water transitions, hydrogen transitions, and rates dictionary.
    """
    data = pd.read_csv(filename, skiprows=17, sep='\s+', header=None)
    temperature = data.iloc[0, 4:].values.astype(float)
    col_names = ['IS_H2O', 'FS_H2O', 'IS_H2', 'FS_H2'] + [f'T_{temp}' for temp in temperature]
    data.columns = col_names
    transition_water    = data.iloc[1:,0:2].values.astype(int) # IS_H2O and FS_H2O
    transition_hydrogen = data.iloc[1:,2:4].values.astype(int) # IS_H2  and FS_H2
    rates_dict = {temp: data.iloc[1:,4+i].values.astype(float) for i, temp in enumerate(temperature)} # Create a dictionary of state-to-state rates for each temperature. Key - temp and Value - rates.
    return temperature, transition_water, transition_hydrogen, rates_dict

# Read the energy levels of hydrogen from the file.
def read_energy_levels_hydrogen(filename):
    """
    Read energy levels of hydrogen from a file.
    
    Args:
    filename (str): Path to the input file.
    
    Returns:
    dict: Energy levels for each J value.
    """
    data = pd.read_csv(filename, skiprows=10, sep='\s+', header=0)
    ground_state = data.groupby('J').first().reset_index()   # Get the ground state energy levels for each J 
    energy_levels = {row['J']: (row['Energy'], row['J'], row['v']) for i, row in ground_state.iterrows()} # Create a dictionary of energy levels for each J. Key - J and Value - (Energy, J, v)
    return energy_levels

# Read the energy levels of water from the file
def read_energy_levels_water(filename):
    """
    Read energy levels of water from a file.
    
    Args:
    filename (str): Path to the input file.
    
    Returns:
    dict: Mapping of level index to J, Ka, Kc values.
    """
    data = pd.read_csv(filename, skiprows=11, sep='\s+', header=0)   
    index_to_jkakc = {row['Level']: (row['J'], row['Ka'], row['Kc']) for _, row in data.iterrows()}
    return index_to_jkakc

# Compute effective rate coefficients
def compute_effective_rate_coefficients(transition_water, transition_hydrogen, rates_dict):
    """
    Compute effective rate coefficients.
    
    Args:
    transition_water (np.array): Water transition states.
    transition_hydrogen (np.array): Hydrogen transition states.
    rates_dict (dict): Dictionary of rates for each temperature.
    
    Returns:
    dict: Effective rates for each temperature and transition.
    """
    effective_rates = {}
    for temp in rates_dict.keys(): # Iterate over each temperature (as in the inp) and compute the effective rate coefficients
        effective_rates[temp] = {}
        unique_transitions = np.unique(np.column_stack((transition_water, transition_hydrogen[:, 0])), axis=0) # Get the unique transitions for water and hydrogen. Basically concatenate the two arrays and get the unique rows
        for trans in unique_transitions: # Iterate over each unique transition and compare it with the transitions in the data and finally sum the rates
            is_h2o, fs_h2o, is_h2  = trans
            mask = (transition_water[:, 0] == is_h2o) & (transition_water[:, 1] == fs_h2o) & (transition_hydrogen[:, 0] == is_h2)
            effective_rates[temp][(is_h2o, fs_h2o, is_h2)] = np.sum(rates_dict[temp][mask]) # Here key is a tuple of the three states and value is the sum of the rates for that transition. Effe rate dict Key - temp, (is_h2o, fs_h2o, is_h2) and Value - sum of rates for that transition
    return effective_rates 

# Compute thermal rate coefficients
def compute_thermal_rates(energy_levels, effective_rates_para, effective_rates_ortho, temperature):
    """
    Compute thermal rate coefficients with fallback logic for missing rates.
    
    Args:
    energy_levels (dict): Energy levels for H2.
    effective_rates_para (dict): Effective rates for para H2.
    effective_rates_ortho (dict): Effective rates for ortho H2.
    temperature (list): List of temperatures.
    
    Returns:
    dict: Thermal rates for each temperature and H2O transition.
    """
    kb = 0.6950386692  # Boltzmann constant in cm^-1/K
    thermal_rates = {}

    j_h2_values = [0, 2, 4, 6, 8] # values of J for p-H2
    j_to_level  = {0:1, 1:1, 2:2, 4:3, 6:4, 8:5} # mapping of J to energy level
    h2_levels   = {j: energy_levels[j] for j in j_h2_values} # energy levels for p-H2

    for temp in temperature:
        thermal_rates[temp] = {} 
        kbt = kb * temp
        zpart = sum((2 * j + 1) * np.exp(-energy / kbt) for energy, j, v in h2_levels.values())        
        unique_water_transitions = set((k[0], k[1]) for k in effective_rates_para[temp].keys()) # Get the unique transitions for water
        for is_h2o_val, fs_h2o_val in unique_water_transitions:
            sum_rates = 0
            for j2a in j_h2_values:
                energy, j, v = h2_levels[j2a]
                weight = (2 * j + 1) * np.exp(-energy / kbt)
                j2a_level = j_to_level[j2a]
                rate = effective_rates_para[temp].get((is_h2o_val, fs_h2o_val, j2a_level))                
                # Fallback logic for missing rates
                if rate is None and j2a == 2:  # If rate is None (zero) for J2a = 2, use J2a = 1 rate from ortho data
                    rate = effective_rates_ortho[temp].get((is_h2o_val, fs_h2o_val, j_to_level[1]))              
                elif rate is None and j2a in [4, 6, 8]: # If rate is None (zero) for J2a = 4, 6, 8, try to use J2a = 2 rate
                    rate = effective_rates_para[temp].get((is_h2o_val, fs_h2o_val, j_to_level[2]))
                    if rate is None:
                        rate = effective_rates_ortho[temp].get((is_h2o_val, fs_h2o_val, j_to_level[1]))

                sum_rates += weight * rate
            thermal_rate = sum_rates / zpart
            thermal_rates[temp][(is_h2o_val, fs_h2o_val)] = thermal_rate

    return thermal_rates

# Compute thermal rate coefficients without fallback
def compute_thermal_rates_regular(energy_levels, effective_rates_para, temperature, h2_sym): 
    """
    Compute thermal rate coefficients for para or ortho H2 without cross-symmetry fallback.
    
    Args:
    energy_levels (dict): Energy levels for H2
    effective_rates_para (dict): Effective rates for para or ortho H2
    temperature (list): List of temperatures
    h2_sym (str): H2 symmetry, either 'para' or 'ortho'
    
    Returns:
    dict: Thermal rates for each temperature and H2O transition
    """
    kb = 0.6950386692  # Boltzmann constant in cm^-1/K
    thermal_rates = {}

    if h2_sym == 'para':
        j_h2_values = [0, 2, 4, 6, 8]            # values of J for p-H2
        j_to_level  = {0:1, 2:2, 4:3, 6:4, 8:5}  # mapping of J to energy level
    elif h2_sym == 'ortho':
        j_h2_values = [1, 3, 5, 7]            # values of J for o-H2
        j_to_level  = {1:1, 3:2, 5:3, 7:4}
    else:
        raise ValueError("h2_sym must be either 'para' or 'ortho'")
    
    h2_levels   = {j: energy_levels[j] for j in j_h2_values} # energy levels 
    for temp in temperature:
        thermal_rates[temp] = {} 
        kbt = kb * temp
        zpart = sum((2 * j + 1) * np.exp(-energy / kbt) for energy, j, v in h2_levels.values())
        unique_water_transitions = set((k[0], k[1]) for k in effective_rates_para[temp].keys())# Get the unique transitions for water
        for is_h2o_val, fs_h2o_val in unique_water_transitions:
            sum_rates = 0
            for j2a in j_h2_values:
                energy, j, v = h2_levels[j2a]
                weight = (2 * j + 1) * np.exp(-energy / kbt)
                j2a_level = j_to_level[j2a]
                rate = effective_rates_para[temp].get((is_h2o_val, fs_h2o_val, j2a_level))        

                if h2_sym == 'para' and rate is None and j2a in [4, 6, 8]:  # Fallback logic for missing rates within the same symmetry
                    rate = effective_rates_para[temp].get((is_h2o_val, fs_h2o_val, j_to_level[2])) # If rate is None (zero) for J2a = 4, 6, 8, try to use J2a = 2 rate 
                elif h2_sym == 'ortho' and rate is None and j2a in [3, 5, 7]:
                    rate = effective_rates_para[temp].get((is_h2o_val, fs_h2o_val, j_to_level[1]))# If rate is None (zero) for J2a = 3, 5, 7, try to use J2a = 1 rate
                if rate is not None:
                    sum_rates += weight * rate

            thermal_rate = sum_rates / zpart
            thermal_rates[temp][(is_h2o_val, fs_h2o_val)] = thermal_rate

    return thermal_rates

# Write thermal rate coefficients to a file
def write_thermal_rates_to_file(thermal_rates, temperature, index_to_jkakc, filename):    
    """
    Write thermal rate coefficients to a file.
    
    Args:
    thermal_rates (dict): Thermal rates for each temperature and transition.
    temperature (list): List of temperatures.
    index_to_jkakc (dict): Mapping of level index to J, Ka, Kc values.
    filename (str): Output file name.
    """
    with open(filename, 'w') as f:
        # Write number of transitions
        num_transitions = len(thermal_rates[temperature[0]])
        f.write(f"! NUMBER OF COLLISIONAL TRANSITIONS\n {num_transitions}\n")  
        f.write(f"! NUMBER OF COLLISION TEMPERATURES\n  {len(temperature)}\n") 
        f.write(f"! THERM. RATE COEFFS UNIT - (cm^3 s^-1)\n")
        f.write('\n')        

        # Write column headers
        temp_headers = ''.join(f"{f'{temp:.2f}':^16s}" for temp in temperature)
        f.write(f"{'IS_H2O':^8s}{'FS_H2O':^9s}{temp_headers}\n")         
        sorted_transitions = sorted(thermal_rates[temperature[0]].keys())  # Sort transitions    

        # Write data
        for i, (is_h2o, fs_h2o) in enumerate(sorted_transitions, 1):
            j_i, ka_i, kc_i = index_to_jkakc.get(is_h2o)
            j_f, ka_f, kc_f = index_to_jkakc.get(fs_h2o)
            rates = [thermal_rates[temp][(is_h2o, fs_h2o)] for temp in temperature] 
            j_i, ka_i, kc_i = int(j_i), int(ka_i), int(kc_i)
            j_f, ka_f, kc_f = int(j_f), int(ka_f), int(kc_f) 

            f.write(f"{is_h2o:6d}{fs_h2o:6d}" + ''.join(f"{rate:16.7E}" for rate in rates) + '\n')            
            # f.write(f"{j_i:3d}{ka_i:3d}{kc_i:3d}{j_f:3d}{ka_f:3d}{kc_f:3d}" +  # uncomment if you want to write J Ka Kc values for each transition
            #         ''.join(f"{rate:16.7E}" for rate in rates) + '\n')


# write effective rate coefficients to a file        
def write_effective_rates_to_file(effective_rates, temperature, filename):
    """
    Write effective rate coefficients to a file.
    
    Args:
    effective_rates (dict): Effective rates for each temperature and transition.
    temperature (list): List of temperatures.
    filename (str): Output file name.
    """
    with open(filename, 'w') as f:
        num_transitions = len(effective_rates[temperature[0]])
        f.write(f"! NUMBER OF COLLISIONAL TRANSITIONS\n {num_transitions}\n")
        f.write(f"! NUMBER OF COLLISION TEMPERATURES\n  {len(temperature)}\n")
        f.write(f"! EFF. RATE COEFFS UNIT - (cm^3 s^-1)\n\n")

        # Write column headers
        temp_headers = ''.join(f"{f'{temp:.2f}':^16s}" for temp in temperature)
        f.write(f"{'IS_H2O':^8s}{'FS_H2O':^9s}{'IS_H2':^8s}{temp_headers}\n") 

        # Sort transitions
        sorted_transitions = sorted(effective_rates[temperature[0]].keys())
        # Write data
        for is_h2o, fs_h2o, is_h2 in sorted_transitions:
            rates = [effective_rates[temp][(is_h2o, fs_h2o, is_h2)] for temp in temperature]
            f.write(f"{is_h2o:7d}{fs_h2o:7d}{is_h2:7d}" + ''.join(f"{rate:16.7E}" for rate in rates) + '\n')

# main program
def main():
    # Input state-to-state rate-coeff file names
    filename_para_para = 'Dubernet_Data/H2O-para_H2-para_statetostate_ID148_v1.kij'
    filename_para_ortho = 'Dubernet_Data/H2O-para_H2-ortho_statetostate_ID149_v1.kij'
    filename_ortho_para = 'Dubernet_Data/H2O-ortho_H2-para_statetostate_ID94_v1.kij' 
    filename_ortho_ortho = 'Dubernet_Data/H2O-ortho_H2-ortho_statetostate_ID147_v1.kij'
    # Energy levels file names
    para_h2_energy_levels_file = 'Dubernet_Data/H2-para_IDT4_vt1_IDC148_vc1.lev'
    ortho_h2_energy_levels_file = 'Dubernet_Data/H2-ortho_IDT16_vt1_IDC149_vc1.lev'
    para_h2o_energy_levels_file = 'Dubernet_Data/H2O-para_IDT14_vt1_IDC149_vc1.lev'
    ortho_h2o_energy_levels_file = 'Dubernet_Data/H2O-ortho_IDT1_vt1_IDC94_vc1.lev'

    # Indicate output file names of effective rates
    output_file_eff_para_para =  'output-effective-rates/effective_rates_para_h2o_para_h2.dat'
    output_file_eff_para_ortho = 'output-effective-rates/effective_rates_para_h2o_ortho_h2.dat' 
    output_file_eff_ortho_para = 'output-effective-rates/effective_rates_ortho_h2o_para_h2.dat'
    output_file_eff_ortho_ortho= 'output-effective-rates/effective_rates_ortho_h2o_ortho_h2.dat'   

    # Indicate output file names of thermal rates
    output_file_thermal_para_h2o_para_h2 = 'output-thermal-rates/thermal_rates_output_para_h2o_para_h2.dat' 
    output_file_thermal_para_h2o_ortho_h2 = 'output-thermal-rates/thermal_rates_output_para_h2o_ortho_h2.dat' 
    output_file_thermal_ortho_h2o_para_h2 = 'output-thermal-rates/thermal_rates_output_ortho_h2o_para_h2.dat' 
    output_file_thermal_ortho_h2o_ortho_h2 = 'output-thermal-rates/thermal_rates_output_ortho_h2o_ortho_h2.dat'  
    
    # Read state-to-state data for para and ortho H2
    temperature, transition_water_pp, transition_hydrogen_pp, rates_dict_para_para      = read_state_to_state_data(filename_para_para)
    temperature, transition_water_po, transition_hydrogen_po, rates_dict_para_ortho     = read_state_to_state_data(filename_para_ortho)
    temperature, transition_water_op, transition_hydrogen_op, rates_dict_ortho_para     = read_state_to_state_data(filename_ortho_para)
    temperature, transition_water_oo, transition_hydrogen_oo, rates_dict_ortho_ortho    = read_state_to_state_data(filename_ortho_ortho) 

    # Compute effective rates for both para and ortho H2
    effective_rates_para_para   = compute_effective_rate_coefficients(transition_water_pp, transition_hydrogen_pp,  rates_dict_para_para)
    effective_rates_para_ortho  = compute_effective_rate_coefficients(transition_water_po, transition_hydrogen_po, rates_dict_para_ortho)
    effective_rates_ortho_para  = compute_effective_rate_coefficients(transition_water_op, transition_hydrogen_op,  rates_dict_ortho_para)
    effective_rates_ortho_ortho = compute_effective_rate_coefficients(transition_water_oo, transition_hydrogen_oo, rates_dict_ortho_ortho)
    
    # Write effective rates to files
    write_effective_rates_to_file(effective_rates_para_para,  temperature,  output_file_eff_para_para)
    write_effective_rates_to_file(effective_rates_para_ortho, temperature,  output_file_eff_para_ortho) 
    write_effective_rates_to_file(effective_rates_ortho_para, temperature,  output_file_eff_ortho_para) 
    write_effective_rates_to_file(effective_rates_ortho_ortho, temperature, output_file_eff_ortho_ortho) 

    # Read energy levels and compute thermal rates
    para_h2_energy_levels =  read_energy_levels_hydrogen(para_h2_energy_levels_file)
    ortho_h2_energy_levels = read_energy_levels_hydrogen(ortho_h2_energy_levels_file)

    thermal_rates_para_para  = compute_thermal_rates(para_h2_energy_levels, effective_rates_para_para, effective_rates_para_ortho, temperature)

    # compute thermal rates without fall back logic to fill missing rates
    thermal_rates_para_ortho  = compute_thermal_rates_regular(ortho_h2_energy_levels, effective_rates_para_ortho, temperature,  h2_sym='ortho') 
    thermal_rates_ortho_para  = compute_thermal_rates_regular(para_h2_energy_levels, effective_rates_ortho_para, temperature,  h2_sym='para')
    thermal_rates_ortho_ortho = compute_thermal_rates_regular(ortho_h2_energy_levels, effective_rates_ortho_ortho, temperature, h2_sym='ortho')
    
    # Read H2O energy levels for J Ka Kc conversion
    para_h2o_index_to_jkakc = read_energy_levels_water(para_h2o_energy_levels_file)
    ortho_h2o_index_to_jkakc = read_energy_levels_water(ortho_h2o_energy_levels_file)
    
    # Write thermal rates to files
    write_thermal_rates_to_file(thermal_rates_para_para, temperature,  para_h2o_index_to_jkakc, output_file_thermal_para_h2o_para_h2)
    write_thermal_rates_to_file(thermal_rates_para_ortho, temperature, para_h2o_index_to_jkakc, output_file_thermal_para_h2o_ortho_h2)
    write_thermal_rates_to_file(thermal_rates_ortho_para, temperature,  ortho_h2o_index_to_jkakc, output_file_thermal_ortho_h2o_para_h2)
    write_thermal_rates_to_file(thermal_rates_ortho_ortho, temperature,  ortho_h2o_index_to_jkakc, output_file_thermal_ortho_h2o_ortho_h2)   

if __name__ == '__main__': 
    main() 