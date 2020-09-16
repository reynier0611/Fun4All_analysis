# Fun4All analysis

This repository can be used to quickly analyze output root files from Fun4All in general, but specifically from [here](https://github.com/eic/g4lblvtx "EIC prototype All-Silicon tracker repository"). This repository is structured in two separate directories: 1) processing, and 2) analysis.

## processing

The first step is done in this repository. Codes here loop over root files expected to be found in a directory called *data*, which should be created inside the top directory. These codes produce secondary root files with extracted momentum, angular, vertex, ... resolutions.

## analysis

The second step is to write codes in this directory to analyze and plot the secondary root files produced in the *processing* step.
