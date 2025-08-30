# solarsystem-mercuryprecession-simulatio

Python simulations of Solar System dynamics and Mercury’s relativistic orbital precession  

---

## Overview  

This project demonstrates two astrophysical simulations implemented from scratch in Python:  

1. **Solar System Simulation** – a 3D N-body model of planetary orbits (including Earth–Moon system).  
2. **Mercury Precession Simulation** – models the perihelion precession of Mercury’s orbit using relativistic corrections to Newtonian gravity.  

The simulations use **numerical integration** and **vectorized physics calculations** with NumPy, and produce animated and plotted outputs with Matplotlib.  

This project is both a showcase of **orbital mechanics** and an example of applied **scientific simulation modeling**.  

---

## Features  

- N-body gravitational dynamics with Newtonian gravity  
- 3D animated visualization of planetary orbits  
- Earth–Moon system modeled with orbital inclination  
- Relativistic correction for Mercury’s perihelion precession  
- Angular perihelion shift analysis in arcseconds per orbit  
- Results consistent with observed ~6 arcseconds/year  

---

## Technology Stack  

- **Python**: NumPy, Matplotlib, mpl_toolkits.mplot3d  
- **Scientific Computing**: Numerical integration, vector math  
- **Visualization**: Matplotlib 3D & animation  

---

## Installation  

Clone the repository:  

```bash
git clone https://github.com/BenMontalbano/solar-system-simulator.git
cd solar-system-simulator
pip install numpy matplotlib
