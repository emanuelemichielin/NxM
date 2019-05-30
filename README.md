This software implements an algorithm to reconstruct the phonon signal from detectors with the SuperCDMS SNOLAB sensor geometry, based on the NxM filter. It includes programs to:

1. Calculate the cross-correlation function of the noise.
2. Calculate templates from data.
3. Apply the NxM filter.

The only requirements are numpy and pyROOT (available on galba.stanford.edu). If working on galba.stanford.edu, the environment can be set up by sourcing the file ``setup_galba.sh``,

``$ source setup_galba.sh``

A simple Monte Carlo example is available under the directory ``example_simple/``, that can be executed as follows:

``$ cd example_simple``
``$ python run.py``

In general, the user will be able to use this example for his own purposes by doing the following modifications:

1. Adapt ``example_simple/DataReader.py`` in order to read the data files of interest.
2. Re-define the radial and azimuthal partitions in ``core/calc_r.py`` and ``core/calc_theta.py``.
3. Re-define the segmentation of the partition space in ``example_simple/bins_part.py``.

A more advanced example is available under the directory ``example_best_mctruth/``, that can be executed as follows:

``$ cd example_best_mctruth``
``$ python run.py``

Some important classes are:

1. ``core/NoisePSDGenerator.py``, that calculates the cross-correlation function of the noise.
2. ``core/TemplateGeneratorNxM.py``, that calculates the templates from data.
3. ``core/OFManagerNxM.py``, that handles the pre-filter step.
4. ``core/OptimalFilterNxM.py``, that applies the NxM filter.

For additional information please see the comments in the code, and the note http://titus.stanford.edu/cdms_restricted/elias/2018_12_06-nxm_filter/index.html.
