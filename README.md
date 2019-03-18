# Beam-Analysis
Automates the process of taking in beam visibilities and creating plots, fits, and maps

INPUTS:

What you have to prepare before using this code:
1)Prepare a directory to use as the current working directory that
  contains the data on all of the baselines. The data should be split
  by imaginary and real values, with naming conventions like:
"I_246_3srcNP_20180101234415_20180102004415"
2)Call main, giving it the "identifier", which in the previous example
would be the 3srcNP, the timestamps in a list (just the unique ones),
the frequencies that data was taken for (in a list), and the seconds
that should be removed either because the data wasn't for the beam in
question or calibrations were occurring during those times.

OUTPUTS:

1)an errors.txt that contains any errors that may have occurred during
this process. (that I expected that I may have to catch)
