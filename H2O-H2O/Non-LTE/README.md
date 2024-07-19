<h1 align="center">H<sub>2</sub>O + H<sub>2</sub>O Collisional Rates Coefficients</h1>

## Objective:

This program computes rate coefficients for rotational state-to-state transitions in H<sub>2</sub>O + H<sub>2</sub>O collisions as a function of rotational and kinetic temperatures of water molecules, *T*<sub>rot</sub> and *T*<sub>kin</sub>, and for any value of the ortho/para ratio of water. The code is primarily developed to be used for astronomical modeling of cometary comae, atmospheres of icy planets, water rich exoplanets and other astrophysical environments where physical conditions deviate from thermodynamic equilibrium and where the H<sub>2</sub>O + H<sub>2</sub>O collisions are important for energy transfer. The **temperatures (both *T*<sub>rot</sub> & *T*<sub>kin</sub>)** in the units of **Kelvin** are used as input and the rotational state-to-state transition **rate coefficients (*k*)** in the units of **cm<sup>3</sup>s<sup>-1</sup>** are given as output. For the target H<sub>2</sub>O molecules **231 transitions** between lower energy para-states and **210 transitions** between lower energy ortho-states are included. For the projectile H<sub>2</sub>O the value of ortho/para ratio can be **from 0 to 1**. The range of temperatures is **5 ≤ *T* ≤ 1000 K**.

## Installation:

This code is written in Fortran language. The only requirement to compile and use this code is to install a Fortran compiler, such as gfortran or ifort. Here are the steps to compile this code and use for the astronomical modeling.

1. First, one needs to download this project by from the [GitHub website](https://github.com/bikramaditya-mandal/Water_Rate_Coefficients.git). Alternatively, one can clone this project using CLI and the following commands.

```sh
   git clone https://github.com/bikramaditya-mandal/Water_Rate_Coefficients.git
   cd Water_Rate_Coefficients
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

## Input:

In the example provided the input values of tempetatures (both *T*<sub>rot</sub> & *T*<sub>kin</sub>) are specified in the file [**Generate_Rates.f90**](Generate_Rates.f90) and are set to *T*<sub>rot</sub> = 250 K and *T*<sub>kin</sub> = 350 K. They can be easily changed to your desired input values of  rotational and kinetic temperatures. Here are the steps to insert the desired temperatures.

1. First, open the file [**Generate_Rates.f90**](Generate_Rates.f90) using any editor.
2. Then, scroll down a bit to locate the line

```sh
    Temp_rot = 250.d0
```

and set the value of this variable to the desired rotational temperature.

**Note:** One can modify it further to run a loop over temperature through the desired range. For example:

```sh
    do i = 1, 10
        Temp_rot = i * 10.0d0
        .
        .
        .
    end do
```

In this example, the code will run calculations for 10 values of rotational temperature: 10, 20, ..., 100 K.

3. Following a similar approach the value of kinetic temperatures can also be modified in the following line:

```sh
    Temp_kin = 350.d0
```

**Note:** One can also add a loop over kinetic temperatures. For example:

```sh
    do i = 1, 3
        Temp_rot = i * 10.0d0
        .
        .
        .

        do j = 1, 5
            Temp_kin = j * 50.0d0
            .
            .
            .
        end do    ! Ending the loop over the kinetic temperature
    end do    ! Ending the loop over the rotational temperature
```

In this example, the code will compute rate coefficients at 3 values of rotational temperature: 10, 20, and 30 K, and, in each case 5 values of kinetic temperature will be considered: 50, 100, 150, 200, and 250 K. 

**Note:** Both *T*<sub>rot</sub> and *T*<sub>kin</sub> can be varied through the range 5 ≤ *T* ≤ 1000 K.

## Output:

After the temperatures are set to the desired values, one would need to recompile the code by using the following commands:

```sh
    make clean
    make
```

To obtain the rate coefficients, one would need to simply run the executable by using the following commands:

```sh
    ./compute_rate.exe
```

Upon execution, the code prints output for a total of 441 rotational state-to-state transitions between first 22 para-H<sub>2</sub>O and 21 ortho-H<sub>2</sub>O states on the screen.

**Note:** For advanced users, the code can be modified based on the need to print into a file for multiple values of temperatures. Alternatively, it can be modified to be used as function which returns the rate coefficients for given input values of temperatures. However, this is recommended for advanced users, and should be done correctly to produce correct rate coefficients.

## Citing this work:

For more details and to cite this work, please refer to:
1. Bikramaditya Mandal et al, Rotational state-to-state transition rate coefficients for H<sub>2</sub>O + H<sub>2</sub>O collisions at nonequilibrium conditions, Astronomy & Astrophysics, 2024.
2. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 671, A51.
3. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 678, A51.
 



