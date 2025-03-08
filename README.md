# Nautilus Submarine Propeller Shaft Analysis - README.md
## Description

This MATLAB project analyzes the structural integrity of a propeller shaft for a Nautilus-class nuclear submarine. The project calculates shear force, bending moment, torque, and angle of twist along the shaft's axis. It also generates design charts for factor of safety, shaft twist ratio, and mass based on material properties for various shaft diameters using a modified Goodman analysis.

## Key Calculations

The project calculates the following key parameters:

*   Shear Force (V)
*   Bending Moment (M)
*   Torque (T)
*   Angle of Twist (Phi)
*   Axial Stress (SigmaX)
*   Mean Stress (SigmaM)
*   Alternating Stress (SigmaAR)
*   Factor of Safety (Xsm)

## Variables

*   `Pprop`: Propeller power requirement (W)
*   `Pgen`: Generator Power Requirement (W)
*   `D1`: Shaft Diameter D1 (m)
*   `D2`: Shaft Diameter D2 (m)
*   `E`: Young's Modulus (N/m2)
*   `pois`: Poisson's Ratio (-)
*   `dens`: Material density (kg/m3)
*   `SigmaU`: Ultimate stress (N/m2)
*   `Nf`: Desired Shaft Life (Cycles)
*   `L1`, `L2`, `L3`, `L4`, `L5`, `L6`: Axial Locations (m)
*   `Rg`: Gear Radius (m)
*   `v`: Submarine Linear Speed (km/hr)
*   `w`: Propeller Shaft Rotational Speed (Rad/s)
*   `Nf1`: Cycles for 0.9SigmaU (S-N Curve)
*   `Nf2`: Cycles for 0.5SigmaU (S-N Curve)

## Equations

The script implements the following equations:

Torque at the propeller
* Tout = Pprop/w

Torque at the gear
* Tgear = Pgen/w

Shear Force
* q(z) = R1z<x-L1>^-1 + R2z<x-(L1+L2+L3)>^-1 + Fz<x-(L1+L2+L3+L4)>^-1
* V(z) = R1z<x-L1>^0 + R2z<x-(L1+L2+L3)>^0 + Fz<x-(L1+L2+L3+L4)>^0
* M(z) = R1z<x-L1>^1 + R2z<x-(L1+L2+L3)>^1 + Fz<x-(L1+L2+L3+L4)>^1

## MATLAB Code Overview

The MATLAB code performs the following main steps:

1.  **Initialization:** Sets initial conditions, material properties, and submarine parameters.
2.  **Calculations:** Calculates shear forces, bending moments, torque, and twist angle along the shaft.
3.  **Stress Analysis:** Performs stress analysis based on the calculated loads and the material properties.
4.  **Goodman Analysis:** Computes factor of safety using a modified Goodman criterion.
5.  **Design Charts:** Generates plots showing the factor of safety, shaft twist ratio, and mass as a function of shaft diameter.
## Dependencies
*   MATLAB

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
