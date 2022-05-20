# semicollisional_electron_tails_modes
Software for finding linear eigenmodes in a reduced, semi-collisional, model for electron-driven microinstabilities in toroidal magnetic confinement fusion devices.


To install, download the program files and install python version 3.6.9.

Create a virtual environment and install the required packages from the requirements.txt file. For instructions, see, for example, https://docs.python.org/3/tutorial/venv.html

To run the eigensolver in the main directory use 

$ python3 sct_run.py

To run the postprocessing diagnostics only use

$ python3 run_post_processing.py

To test the transport coefficients use

$  python3 trapped_fraction_test.py

This program currently uses the geometry output from the gyrokinetic code GS2 (https://bitbucket.org/gyrokinetics/gs2/src/master/) to specify the equilibrium magnetic field.
An example GS2 input and output file is provided. In addition, an option to specify the coefficients in the equations manually is provided in the source code. 

Documentation on the equations solved, and their derivation, is pending.
