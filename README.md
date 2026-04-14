# vasco
VASCO: VLBI and SMILE source based CASA Optimizations

This is a helper tool for CASA based optimizations to use rPICARD and CASA for VLBI in general and in particular for the SMILE sources.


> UNDER DEVELOPMENT : Version 0.3.0 releasing soon.

# Installation

```
git clone --recurse-submodules https://github.com/avikhagol/vasco.git
```

# Usage

```
vasco pipeline run [OPTIONS] [FITSFILENAMES]... 
```

```
vasco listobs [FITSFILENAMES]... 
```


# Development
After you clone the repo
```bash

$ cd vasco
$ pip install -e .[dev]
```

# Attribution

When using VASCO, please add a link to this repository in a footnote.

# Acknowledgement

"VASCO was developed within the "Search for Milli-Lenses" (SMILE) project. SMILE has received funding from the European Research Council (ERC) under the HORIZON ERC Grants 2021 programme (grant agreement No. 101040021).