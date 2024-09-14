<h1 align="center">H<sub>2</sub>O + H<sub>2</sub>O Collisional Rates Coefficients</h1>

## Objective:

This program computes rate coefficients for state-to-state rotational transitions in H<sub>2</sub>O + H<sub>2</sub>O system for a Local Thermodynamic Equilibrium (LTE) scenario where Rotational and Kinetic temperature (T<sub>rot</sub> & T<sub>kin</sub>) is equal to each other (i.e., T<sub>rot</sub> = T<sub>kin</sub>). The code is primarily developed to be used for astronomical modeling in atmospheres of comets, and other astrophysical environments where H<sub>2</sub>O is discovered. Only one **temperature (representing both T<sub>rot</sub> & T<sub>kin</sub>)** is used as input in the unit of **kelvin** and the rotational state-to-state transitions **rate coefficients (*k*)** is computed as ***cm<sup>3</sup>s<sup>-1</sup>***.

## Installation:

This code is written in Fortran language. The only requirement to compile and use this code is to install a Fortran compiler, such as gfortran or ifort. Here are the steps to compile this code and use for the astronomical modeling.

1. First, one needs to download this project from the [GitHub website](https://github.com/bikramaditya-mandal/Rate-Coefficients.git). Alternatively, one can clone this project using CLI and the following commands.

```sh
   git clone https://github.com/bikramaditya-mandal/Rate-Coefficients.git
   cd Rate-Coefficients/H2O-H2O/LTE/
```

2. Then, one needs to edit the first line of the Makefile to incorporate the appropriate installed compiler by modifying 

```sh
    FC=gfortran
```

and  replace **gfortran** with the choice of compiler.

3. Then one needs to clean the directory to remove old and unnecessary module and object files using 

```sh
    make clean
```

4. Finally, use "make" to compile and create the executable.
```sh
    make
```

## Modify temperatures:

There is only one file which needs to be modified to get rate coefficients for the desired Temperature (since T<sub>rot</sub> = T<sub>kin</sub>). The file name is [**Generate_Rates.f90**](Generate_Rates.f90). In this example file, both rotational and kinetic temperatures are set at ***Temp = 300 k***.

Here are the steps to insert desired temperatures.

1. First, open the file [**Generate_Rates.f90**](Generate_Rates.f90) using any editor.
2. Then, scroll down a bit to locate the line

```sh
    Temp = 300.d0
```

and set the value of the desired **Rotational Temperature** to this variable.

**Note:** One can modify it further to run a loop over ranges of temperatures. For example:

```sh
    do i = 1, 10
        Temp = i * 10.0d0
        .
        .
        .
    end do
```

In this example, the code will compute 10 values of **Temperatures** of **10, 20, ..., 100 *k***.


## Output:

After the temperatures are set to the desired values, one would need to recompile the code by using the following commands:

```sh
    make clean
    make
```

To get the rate coefficients, one would need to simply execute the executable by using the following commands:

```sh
    ./compute_rate.exe
```

Upon execution, the code prints output for a total of 441 rotational state-to-state transitions between first 22 *para*-H<sub>2</sub>O and 21 *ortho*-H<sub>2</sub>O states on screen.

**Note:** For advanced users, the code can be modified based on the need to print into a file for multiple values of temperatures. Alternatively, it can be modified to be used as function which returns the rate coefficients for given values of temperatures. However, this is recommended for advanced users, and should be done correctly to produce correct rate coefficients.

## Citing this work:

For more details and to cite this work, please refer to:
1. Bikramaditya Mandal et al, 2024, Astronomy & Astrophysics, 688, A208.
2. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 671, A51.
3. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 678, A51.

**Note:** This code works for **Rotational Temperature = Kinetic Temperature** in the range **5&le;*Temp*&le;1000 *k***. 


