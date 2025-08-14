# **H2O-H2 Collision Cross-section Data**
Collision cross-sections for H2O-H2 molecular systems at ten collision energies.

### **Directory Structure**
The data is organized by molecular symmetry combinations, with each ZIP file containing collision cross-section data across all 10 collision energies:
- **`para-H2O_para-H2.zip`** - Para-H2O + Para-H2 collisions
- **`para-H2O_ortho-H2.zip`** - Para-H2O + Ortho-H2 collisions
- **`ortho-H2O_para-H2.zip`** - Ortho-H2O + Para-H2 collisions
- **`ortho-H2O_ortho-H2.zip`** - Ortho-H2O + Ortho-H2 collisions

  
### **Data Format**
- **j_i, ka_i, kc_i → j_f, ka_f, kc_f**: H2O quantum numbers (quenching direction)
- **H2_i(j) → H2_f(j)**: H2 rotational states
- **sigmaU_qu, sigmaE_qu**: Quenching cross-sections 
- **sigmaU_ex, sigmaE_ex**: Excitation cross-sections 
- **E_coll_qu, E_coll_ex**: Respective collision energies
- 
**Note:** Quantum numbers correspond to quenching transitions. Excitation cross-sections are for the reverse process.
