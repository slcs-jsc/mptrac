# Quickstart

A simple example is provided, illustrating how to simulate the dispersion of volcanic ash from the eruption of the Puyehue-Cordón Caulle volcano, Chile, in June 2011.

## Example of an MPTRAC simulation

The example can be found in the project directory. The project directory can also be used to store results of other simulation and experiments with MPTRAC. The simulation is controlled by a shell script:

    cd mptrac/projects/example
    ./run.sh

Please see the script `run.sh` on how to invoke MPTRAC programs such as `atm_init` and `atm_split` to initialize trajectory seeds and `trac` to calculate the trajectories.

The script generates a number of plots of the simulation output at different time steps after the eruption by means of the `gnuplot` graphing tool. These plots should look similar to the output already provided in the repository.

This is an example showing the particle output on 6th to 8th of June 2011:
![MPTRAC particle output for Puyehue simulation](img/puyehue.png)

## Further information

More detailed information for new users and developers is provided in the [MPTRAC user manual](https://slcs-jsc.github.io/mptrac/) and collected in the [GitHub wiki](https://github.com/slcs-jsc/mptrac/wiki).

A detailed description of the MPTRAC model is provided in this paper:

Hoffmann, L., Baumeister, P. F., Cai, Z., Clemens, J., Griessbach, S., Günther, G., Heng, Y., Liu, M., Haghighi Mood, K., Stein, O., Thomas, N., Vogel, B., Wu, X., and Zou, L.: Massive-Parallel Trajectory Calculations version 2.2 (MPTRAC-2.2): Lagrangian transport simulations on graphics processing units (GPUs), Geosci. Model Dev., 15, 2731–2762, https://doi.org/10.5194/gmd-15-2731-2022, 2022.

We are interested in sharing MPTRAC for operational and research applications. Please do not hesitate to contact us, if you have any further questions or need support.
