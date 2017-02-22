Expressive Dynamic Modeling for Violin
--------------------------------------

This file describes the necessary steps to reproduce the experiment of expressive dynamic modeling.

### Dependencies Installation

The resource folder contains all additional software required for running the code, apart from Mathworks MATLAB itself.

#### pYin Vamp Plugin

To install the plugin, unzip the .gz file and copy both pYin.dylib and pYin.cat to the /Library/Audio/Plug-Ins/VAMP folder creating it if it does not exist. (You may need to be a machine administrator for this.)

### Credits

#### pYin Vamp Plugin

The system uses the "notes" algorithm implemented in the pYin Vamp plugin, developed by Matthias Mauch @ QMUL and described in the [Tony paper].

#### Sonic Annotator

The Sonic Annotator is a command-line tool for executing Vamp plugins developed and maintained by QMUL.

#### Midi Toolbox

The MATLAB code uses functions of the Midi Toolbox for reading/writing MIDI files and also its implementation of Emilios Camboroupoulos Local Boundary Detection Model as described in the [LBDM paper].