Dependencies are ROOT, MGDO, and jsoncpp.  Compile with:
./pmake.py

Executables are built to ./bin.  Running them with the -h flag will provide more information on the available options.

To analyze waveforms from a build data file:
bin/process -i input.root -o output.root -j config.json -b 16 -c 0 -w 1000
Where:
  -c specifies the channel index
  -b specifies the number of bits of the ADC (just used for binning histograms
  -j specifies the name of a json configuration file (defaults are assumed if this is not specified)
  -o specifies the output file name
  -i specifies the input file name
  -w sets the interval for writing example waveforms with transforms and time points to the output file
The output file will contain a few diagnostic histograms in addition to a ROOT tree with the values for each event.

With one or a number of outputs from the above, a number of calibration constants are determined from:
bin/calibrate -i input_0.root -i input_1.root -o output.root -c 0 -s Th
where:
  -s specifies the source type (currently only Th or Co)
  -c specifies the channel index
  -o specifies the output filename
  -i specifies the input filename with multiple instances joined together
This will print to the std out a set of calibration constants computed.  If this is the first time this detector has been analyzed, then a new configuration file for the first step above should be created.  The PZ decay time should then be set accordingly there and fed back into the first step of the analysis.

The calibration constants are also written by default to ./calibrate.json (the filename can be set with the -J flag).  In the calibration step, to read in these constants, specify the filename with the -j flag, and the constants will be applied without recomputing the parameters.

As an example of plotting things in the output tree, one can run:
root -l 'plot_calibrated.C("calibrated_file.root")'
This will produce a few histograms as examples.  Note that the rise time correction is a work in progress, and the AvsE calibration also needs some attention.

pulse_shape.cc is a work in progress to extract ADC overshoot parameters and other pulse shape effects.