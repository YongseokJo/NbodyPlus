## This is an ABYSS code simulating self-consistent gravity for stars and black holes.
# How to run this code

> `e.g., abyss.exe -c config.txt`  
<pre>
-c|--config [name of a configuration (parameter) file]
</pre>

>`structure of a configuration file`
<pre>
Filename        = [name of initial condition]
eta             = [constant for the timestepper]
FixNumNeighbor  = [number of neighbors (N^(1/2) is recommended)]
InitialRadius   = [initial radius of neighbor (in pc)]
StopTime        = [duration of simulations (in yr)]
dtOutput        = [frequency of output (in yr)]
OutputDirectory = [name of output directory]
</pre>

# How to compile this code
CUDA is mendatory as of now.

## On-going Projects
- GPU optimization by Minyong Jung
- Few-body dynamics & Post-Newtonian treatment by Eunwoo Chung
- CPU parallelization by Yongseok Jo
- Integration with a hydrodynamic code (Enzo-Abyss) by Yongseok Jo

## Future Projects
- Tidal disruption for black hole
- Gravitational wave
- X-ray Binary (XRB)
- Cosmic ray
- Developement of engines for acceleration computation

For these projects, if interested, please feel free to contact Yongseok Jo (g.kerex@gmail.com).
  

## Contributions ###
* Yongseok Jo     @ Columbia University
* Eunwoo Chung    @ Seoul National Univeristy
* Minyong Jung    @ Seoul National University
* Seoyoung Kim    @ University of Wisconsin-Madison


### Thanks to
* Ji-hoon Kim    @ Seoul National University
* Greg L. Bryan  @ Columbia University

## Contact

Yongseok Jo @ [yj2812@columbia.edu](mailto:yj2812@columbia.edu)

Enzo-Abyss repository : [https://github.com/YongseokJo/enzo-nbody](https://github.com/YongseokJo/enzo-nbody)
