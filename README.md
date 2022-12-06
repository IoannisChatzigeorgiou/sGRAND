This repository contains the MATLAB&reg; code that was used to run simulations and obtain results presented in the paper:
"Symbol-Level GRAND for High-Order Modulation over Block Fading Channels" by Ioannis Chatzigeorgiou and Francisco Monteiro. The paper has been accepted for publication in IEEE Communications Letters.

We acknowledge that the methodology of Guessing Random Additive Noise Decoding (GRAND&trade;) has been developed by Ken Duffy, Jiange Li and Muriel Médard. For more details, please visit: <a href="https://www.granddecoder.mit.edu/" target="_blank">www.granddecoder.mit.edu</a>.

## Reproduce the simulation plots using the available simulation data

Set the folder "Simulation_Results" as the Current Folder in MATLAB<sup>&reg;</sup> and run `plot_BLER` in the Command Window to obtain the relevant figure.

## Re-run the simulations to obtain the data needed for the simulation plots

If you do not want to use the readily available simulation results but you prefer to re-run the simulations in order to obtain the MAT files, run `sim_RLC_16QAM_symbol_GRAND_Rayleigh_BLER` in the Command Window.

Please edit the top lines in the file above if you wish to change the system parameters (e.g., code rate, abandonment threshold) or the channel parameters (e.g., the E<sub>b</sub>/N<sub>0</sub>).
