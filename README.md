# crtanalysis
mostly ROOT macros for playing with ntuples output from icaruscode analyzer modules

Setup
  take a look inside of Analyze.C to see what functions are available. these are organized by which trees they access. At the top of the file, you can see a contruction for a TreeManager object - this is where you should set your full file path to the ntuple you wish to analyze

TO RUN:

$ root
$ .L TreeManager.C
$ .L Analyze.C
$ PlotFunctions...
