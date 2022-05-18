# Water Drop and Ripple Simulation
This is a repo for Water Drop and Ripple simulation in C++ with GUI for different functionalities.
## Description
The goal of our project was to find an efficient way to realistically simulate ripples in a water surface. We implemented our water ripples using surface simulation, modeling a water surface as a series of point masses and springs. Water ripples were initiated by exerting displacement forces on certain point masses, and then simulation of wave dynamics was calculated across the surface to propagate the ripples. We cover both an initial sine-wave based approach, as well as a more efficient linear approximation approach for calculating ripple dynamics. Next, we describe a method for simulating water drops hitting the surface via collision detection. Finally, we provide an implementation of the simulation of the waterdrops and reflection/refraction properties of a water that provide a complete rendering of the surface appearance.

## Getting Started

### Installing

```
git clone https://github.com/M1koto/raindrop_simulation/
```

```
mkdir build; cd build
cmake ..
make
```

### Executing program




## Authors

Kenny Liao, Mark Chan, Alan Zhu, Jim Wang
