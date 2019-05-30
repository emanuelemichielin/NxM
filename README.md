This software implements an algorithm to reconstruct the phonon signal from detectors with the SuperCDMS SNOLAB sensor geometry, based on the NxM filter. It includes programs to:

1. Calculate the cross-correlation function of the noise.
2. Calculate templates from data.
3. Apply the NxM filter.

The only requirements are numpy and pyROOT. I reccomend to work on slac cluster, and to get python 3 and root6

"$ source /nfs/slac/g/supercdms/analysis/anaconda3/bin/activate"

To run the code, go to data_analysis and then 

``$ python run.py``

Some important classes are:

1. ``core/NoisePSDGenerator.py``, that calculates the cross-correlation function of the noise.
2. ``core/TemplateGeneratorNxM.py``, that calculates the templates from data.
3. ``core/OFManagerNxM.py``, that handles the pre-filter step.
4. ``core/OptimalFilterNxM.py``, that applies the NxM filter.

The DataReader OpenFile method requires the path and the name of files to process (you can just give the path and it will process everything) and the number of events (if 0 or non-given will process all the events). The first ~15 dumps are randoms, the other are triggered events. Last week I have learned how to get the trigger info from the event, I need to implemented that to be automatic.

A very raw, event-by-event, selection is done to get rid of "bad" events (pile up, muon tails and saturated events). This is absolutely preliminary, and need to be improved. (see autocus of QETpy package for example).

Everything, up to the filter application, seems to work. I would much appreciated if you can start looking at this last step and why it is giving random results (to double check if the templates are appropiate).
