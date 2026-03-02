# Part 2: Zee calibration

Welcome to part 2 of FP2 ATLAS Experiment laboratory. 

To access data, you will have to fix the logical link pointing to ntuples. Please do something like:
```
rm data #(exiting link)
ln -s /path/toyourdirectory/where/ntuple_part23/is/stored data

# Same for fputils
rm fputils #(exiting link)
ln -s ../fputils fputils

```

## Some variables
 ```bash
el_n:            number of electrons per event
el{i}_cl_E:      Uncalibrated energy of the {i}th electron (replace {i} by 1 or 2; the same hereafter)
el{i}_cl_eta:    Eta of the {i}th electron
el{i}_cl_phi:    Phi of the {i}th electron
el{i}_cl_charge: Electric charge of the {i}th electron
 ```

## To plot the distribution 
 ```bash
plt.hist(dataframe["el1_isolation"])
plt.hist(dataframe["el1_isolation"], bins=100)
 ```
