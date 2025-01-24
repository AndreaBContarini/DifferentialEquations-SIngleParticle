# Single Particle Motion Simulation

This repository contains the code and resources for simulating the motion of a single particle under different physical conditions. The project explores the motion governed by second-order differential equations, with and without friction, using numerical integration techniques like Velocity Verlet and Runge-Kutta 2. The study aims to understand the particle's dynamics under varying initial velocities and friction coefficients.

## Project Overview
The motion of a single particle of mass \( m = 1 \) is simulated along a one-dimensional path. The study is divided into three parts:
1. **Motion without friction**: Integration using the Velocity Verlet algorithm.
2. **Motion with friction**: Studying the effect of a friction term using the Runge-Kutta 2 algorithm.
3. **Critical regions analysis**: Identifying regions in the velocity-friction parameter space where the particle stabilizes at specific positions.

## Features
- **Simulation Algorithms**:
  - Velocity Verlet for high-precision integration in frictionless motion.
  - Runge-Kutta 2 for scenarios including friction.
- **Parameter Exploration**:
  - Critical initial velocities (\( v_0 \)) for the particle to reach specific positions.
  - Identification of regions in the \( v_0 \)-\( \gamma \) plane where the particle stabilizes asymptotically.
- **Visualization**:
  - Phase space and trajectory plots for various initial conditions.
  - Region plots to study stabilization under long-term dynamics.

## Key Results
- Demonstrated the second-order accuracy of the Velocity Verlet algorithm.
- Identified critical values of \( v_0 \) and \( \gamma \) where the particle oscillates around or stabilizes at \( x = 1 \).
- Analyzed the impact of integration step size (\( \Delta t \)) on precision.


## Prerequisites
- GCC or a similar C compiler
- Gnuplot for visualization
- Basic knowledge of terminal commands

## Compilation and Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/Single_Particle_Motion_Simulation.git
   cd Single_Particle_Motion_Simulation
