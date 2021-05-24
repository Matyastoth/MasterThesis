# MasterThesis
A collection of methods used during the MSc thesis project, done in Java (NetBeans 8.0.2).

The file containing the scripts can be found in MasterThesis/MasterThesis/src/masterthesis/MasterThesis.java.

During the MSc thesis, my goal was to predict the long-term fouling rate of membranes in an MBR system.
I was working with a year's worth of data from a MBR sewage water cleaning facility; as such, the project
included data culling (by calculating std. deviations), smoothing (5 point Lagrange polynomial interpolation),
recursive fitting (exponential formulas used in this field).

This latter outsourced the hardest part of the work to Matlab by calling it through a system command, and extracting the parameters.
(To be honest, I'm still not sure whether the dev. team intended it to be used such way, but I couldn't be bothered.)
