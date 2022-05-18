# Water Drop and Ripple Simulation
This is a repo for Water Drop and Ripple simulation in C++ with GUI for different functionalities.
## Description
The goal of our project was to find an efficient way to realistically simulate water drops and ripples in a water surface. We implemented our water ripples using surface simulation, modeling a water surface as a series of point masses and springs. Water ripples were initiated by exerting displacement forces on certain point masses, and then simulation of wave dynamics was calculated across the surface to propagate the ripples. We cover both an initial sine-wave based approach, as well as a more efficient linear approximation approach for calculating ripple dynamics. Next, we describe a method for simulating water drops hitting the surface via collision detection. Finally, we provide an implementation of the simulation of the waterdrops and reflection/refraction properties of a water that provide a complete rendering of the surface appearance.

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
```
./clothsim -f ../scene/sphere.json
```

## GUI
Please use the GUI provided to adjust water drop position, simulation speed, damping, and gravity if needed.

<img width="233" alt="Screen Shot 2022-05-18 at 3 02 53 PM" src="https://user-images.githubusercontent.com/60208038/169163340-7eeada1b-5b85-4c29-ab57-3490771e01db.png">

Please use the keyboard commands if you'd like.

<img width="644" alt="Screen Shot 2022-05-18 at 3 07 59 PM" src="https://user-images.githubusercontent.com/60208038/169163949-93e1cb46-0241-4a73-8a9c-b7556b89a03e.png">


## Physical Simulation
Please check the complete physical simulation report at: 
https://docs.google.com/document/d/1HSeB-ANhKU0oZETl8O9rdO01HcSX7cRG4B_7fHjJKYA/edit?fbclid=IwAR0sybMe8KWVDUUvkPaTRShniPv2P2xlUf8WcCP5nY0CuL0wpmg9hKRdpKA#heading=h.j2s3lk3lh9lo

## Demo
### Single Water Drop

https://user-images.githubusercontent.com/60208038/169162801-f0062655-299b-45f6-b0aa-e12e64fd0f8c.mp4

### Multiple Water Drop

https://user-images.githubusercontent.com/60208038/169162981-92f9044d-6f7a-4cb1-8e65-a8ba43bc3f47.mp4

### Complete Demo
https://www.youtube.com/watch?v=-Y__BQRy5Do&ab_channel=JimWang
## Authors

Kenny Liao, Mark Chan, Alan Zhu, Jim Wang

## Credits

CS184 Project4: Clothsim: https://cs184.eecs.berkeley.edu/sp22/docs/proj4

X. Zhang and G. Yang, "Ripple simulation based on mass-spring model and image space computation," 2010 3rd International Congress on Image and Signal Processing, 2010, pp. 553-557, doi: 10.1109/CISP.2010.5647276.
