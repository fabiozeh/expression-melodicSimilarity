Expressive Dynamic Modeling for Violin
--------------------------------------

This repository contains code developed at Universitat Pompeu Fabra for supporting research in music expressivity modelling.

All code available has been developed and tested in Mathworks MATLAB r2016a. No additional installations required. Due to licensing restrictions, not all datasets used can be made available.
To reproduce the experiment of expressive dynamic modeling for violin, refer to the script _main.m_ in the code folder.

When cloning this repository, remember to do so with the "recursive" option (Git 1.6.5 or later) in order to load submodules as well.
```
git clone --recursive https://github.com/fabiozeh/expressiveViolin.git
```

### Credits

Some third party code dependencies are included, namely:

#### Midi Toolbox

The MATLAB code uses functions of the [Midi Toolbox](https://github.com/miditoolbox/1.1) for reading/writing MIDI files and also its implementation of Emilios Camboroupoulos Local Boundary Detection Model as described in the [LBDM paper](http://users.auth.gr/emilios/papers/icmc2001.pdf).
