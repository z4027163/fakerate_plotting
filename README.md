# Ks Fake Rate Measurement Code

## Requirements

This code was created under CMSSW_11_3_4. No other packages needed.

## Running the Fakerate Measurement

The measurement can be decomposed into two different sub tasks

- creating the histograms for fitting
- fitting and creating plots

The first task is done by the `plotFakenew.cc` file and the second is done by the `fitks.cc` file. The `plotFakenew.cc` must first be compiled. To do this, simply run `make -j 8` to run the make file. Then, one can simply run `./Fake` in the terminal to start producing the histograms. You do not need to compile `plotFakenew.cc` seperately, instead it can simply be run by using `root -l -b -q fitks.cc`. Additionally, since the std output to the terminal is quite large, it is recommended to save the output to a log file using a command such as `root -l -b -q fitks.cc >&log_fits&`.

One can additionally run the entire measurement by using the `run_all.sh` script. However, for this to work the inputs and outputs of the two files (`plotFakenew.cc` and `fitks.cc`) must be properly configured. See below for more information.

### Histogramming code (`plotFakenew.cc`)

The inputs to the histogramming code are file lists. The file lists contain directories (ending in `/`) which are seperated by new line characters (no spaces, no commas). The code will then look at every file in each of the directories specified. 

Additionally, in the `main` function, the code requires you to specify output directory under the `path` variable and the binning type under the `binning` variable. Currently the code only supports binning by `pT` and `lxy`. 

Lastly, for each file list which is given, once adds the following 5 lines of code at the end of the file:

```
  plotFakenew c1;
  TChain *tC1 = new TChain("Events");
  c1.list_dir_file(tC1, "NAME OF INPUT FILE LIST HERE");
  cout << "Entry Data " << tC1->GetEntries() << endl;
  c1.loopOverChain(tC1, Form("DESCRIPTIVE_NAME_OF_OUTPUT_FILE_HERE_%s", binning.Data()), "Histproduction", path, binning);
  delete tC1;
```

The output of this code will then be ROOT files, named and located as parametrized above, with histograms.

Note that sometimes, because histogramming code can take some time on large files, it may be benificial to manually created combined datasets from the indvidual ones (for example combininng EGamma and Parking data). One can use the ROOT `hadd` command as show by this example:

```
hadd histo_data_combined_522_2022_ks_trigger_pT.root histo_data_*
```

### Fitting and plotting code (`fitks.cc`)

The fitting and plotting code takes as input the files created in the previous step. All basic changes to the file can be made in the `fitks_loop`

First, make sure that the `path` variable is properly specified and all the files that are loaded are present. (if not, either create histograms for these new files or simply remove them from the plotting and fitting code.). The code then creates many plots for each of the fits and then binned plots for the fakerate. Which binning to run is specified in the `fitks` function at the very end of the file.

Below is an overview of the main plotting and fitting functions.

The `playfit` function goes over each of the bins and creates a fit and a plot for this fit. Note that fits are not garenteed to converge and any final analysis should be inspected by eye for the plots. These functions also output histograms which are binned and contain the final fake rates per bin.

The `overlay#` functions take in the histograms output by playfit and overlay them for a final output. Currently, `overlay1`, `overlay2`, `overlay3` are the only functions supported, so one can only overlay 1,2, or 3 fits. Eventually, a general function should be made. 

In general, the code then takes as input the histograms from `plotFakenew.cc`, executes `playfit` on each one, and then generates various overlays. All plots are outputed in the same directory as input histograms and are saved in pdf and png versions.

## Plotting Individual Variables

There is one last script (`plotVars.cc`), which allows you to look at the distribution over individual variables. Like the histogramming code, this code must also first be compiled using a command such as `make -j 8`. Then, one can run `./Varplot` to generate the plots. 

The inputs and specifications are all configured in the `main` function of the file and are very similar to the histogramming code. 