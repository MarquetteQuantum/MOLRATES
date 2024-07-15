<h1 align="center">H<sub>2</sub>O + H<sub>2</sub>O Collisional Rates Coefficients</h1>

## Objective:

This program computes rate coefficients for state-to-state rotational transitions in H<sub>2</sub>O + H<sub>2</sub>O system as a function of Rotational and Kinetic temperatures (T<sub>rot</sub> & T<sub>kin</sub>). The code is primarily developed to be used for astronomical modeling in atmospheres of comets, and other astrophysical environments where H<sub>2</sub>O is discovered. The **temperatures (both T<sub>rot</sub> & T<sub>kin</sub>)** are used as input in the unit of **kelvin** and the rotational state-to-state transitions **rate coefficients (*k*)** is computed as ***cm<sup>3</sup>s<sup>-1</sup>***.

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

## Modify temperatures:

There is only one file which needs to be modified to get rate coefficients for the desired Temperatures (both T<sub>rot</sub> & T<sub>kin</sub>). The file name is [**Generate_Rates.f90**](Generate_Rates.f90). In this example file, the rotational temperature is set at ***T<sub>rot</sub> = 250 k*** and the kinetic temperature is ***T<sub>kin</sub> = 350 k***.

Here are the steps to insert desired temperatures.

1. First, open the file [**Generate_Rates.f90**](Generate_Rates.f90) using any editor.
2. Then, scroll down a bit to locate the line

```sh
    Temp_rot = 250.d0
```

and set the value of the desired **Rotational Temperature** to this variable.

**Note:** One can modify it further to run a loop over ranges of temperatures. For example:

```sh
    do i = 1, 10
        Temp_rot = i * 10.0d0
        .
        .
        .
    end do
```

In this example, the code will compute 10 values of **Rotational Temperature** of **10, 20, ..., 100 *k***.

3. Following a similar approach, find the line below

```sh
    Temp_kin = 350.d0
```

**Note:** One can also add a loop over the kinetic temperatures as well. For example:

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

In this example, the code will compute 3 values of **Rotational Temperature** of **10, 20, and 30 *k***. For each rotational temperature, the code will compute 5 values of **Kinetic Temperature** of **50, 100, 150, 200, and 250 *k***.

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
1. Bikramaditya Mandal et al, Rotational state-to-state transition rate coefficients for H<sub>2</sub>O + H<sub>2</sub>O collisions at nonequilibrium conditions, Astronomy & Astrophysics, 2024.
2. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 671, A51.
3. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 678, A51.

**Note:** This code works for **Rotational Temperature** in the range **5&le;*T*<sub>rot</sub>&le;1000 *k***, and the **Kinetic Temperature** in the same range of **5&le;*T*<sub>kin</sub>&le;1000 *k***. 



