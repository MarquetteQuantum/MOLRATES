# **H2O-H2 Collision Cross-section Data**

Collision cross-sections for H2O-H2 molecular systems at ten collision energies.

### **Directory Structure**

Each `UXXXX` folder contains collision data for energy U = XXXX (in cm⁻¹):

- `U12000/` - Collision energy 12000 cm⁻¹
- `U1430/` - Collision energy 1430 cm⁻¹
- `U170/` - Collision energy 170 cm⁻¹
- `U20/` - Collision energy 20 cm⁻¹

### **Data Files** (per energy folder)

Each energy directory contains four symmetry combinations:

- **Para-H2O + Para-H2**
- **Para-H2O + Ortho-H2**
- **Ortho-H2O + Para-H2**
- **Ortho-H2O + Ortho-H2**

### **Data Format**

- **j_i, ka_i, kc_i → j_f, ka_f, kc_f**: H2O quantum numbers (quenching direction)
- **H2_i(j) → H2_f(j)**: H2 rotational states
- **sigmaU_qu, sigmaE_qu**: Quenching cross-sections 
- **sigmaU_ex, sigmaE_ex**: Excitation cross-sections 
- **E_coll_qu, E_coll_ex**: Respective collision energies

**Note:** Quantum numbers correspond to quenching transitions. Excitation cross-sections are for the reverse process.
