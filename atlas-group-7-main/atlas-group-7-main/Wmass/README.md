# Part 3: Wmass measurement

Welcome to part 3 of FP2 ATLAS Experiment laboratory. 

To access data, you will have to fix the logical link pointing to ntuples. Please do something like:
```
rm data #(exiting link)
ln -s /path/toyourdirectory/where/ntuple_part23/is/stored data  

# Same for fputils
rm fputils #(exiting link)
ln -s ../fputils fputils

```


