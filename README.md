# vasco
VASCO: VLBI and SMILE source based CASA Optimizations

This is a helper tool for CASA based optimizations to use rPICARD and CASA for VLBI in general and in particular for the SMILE sources.

# Example

```bash
$ vasco -l "SCAN,SOURCE,ANTENNA" test.fits > list.obs
$ vasco path/to/fits --identify-targets
ec071d_1_1.IDI1
science:['J1143+1834']

phase:['J1148+1840']

brightcal:['J0927+3902', 'J1310+3220']

$ vasco --find-refant path/to/fits/
ec071e_1_1.IDI9
ANNAME  STD_TSYS  nRows   Distance
    JB 15.055768   7880 269.581752
    TR 37.115318   6795 205.491414
    SR 42.076019   2717 270.728365
    O6 74.425308   4383 234.552076

$ vasco -S 0851+071 # search for source '0851+071' in the rfc_3d catalog
Category  IVS name  J2000 name hh(RA) mm(RA)     ss(RA) dd(DEC) mm(DEC)     ss(DEC)   D_alp   D_Del    Corr     Obs S_T-  ... X_T-  X_Tot X_u- X_unres U_T-  U_Tot U_u- U_unres K_T-  K_Tot K_u- K_unres Type       Cat
8296        C  0851+071  J0853+0654     08     53  48.189960     +06      54  47.23480      0.11    0.13  -0.065     686       ...       0.466        0.180    -  1.00     -   1.00     -  1.00     -   1.00   Fus  rfc_2023

$ vasco -S 0851+071 -C C_Tot
C_Tot
8296  0.628

```

# USAGE

```bash
usage: vasco [-h] [-l [LIST_OBSERVATION]] [-t] [-r] [-s SPLIT_SOURCE] [--to-B TO_B] [-S SEARCH_SOURCE] [-C COLUMN] [input_file [input_file ...]]

____    ____  ___           _______.  ______   ______   
\   \  /   / /   \         /       | /      | /  __  \  
 \   \/   / /  ^  \       |   (----`|  ,----'|  |  |  | 
  \      / /  /_\  \       \   \    |  |     |  |  |  | 
   \    / /  _____  \  .----)   |   |  `----.|  `--'  | 
    \__/ /__/     \__\ |_______/     \______| \______/  

    VLBI and SMILE source based CASA Optimizations (VASCO).                                                        
    
    

positional arguments:
  input_file            Give the input file path.

optional arguments:
  -h, --help            show this help message and exit

operations:
  
                               use operations based on file type e.g., .FITS .MS

  -l [LIST_OBSERVATION], --list-observation [LIST_OBSERVATION]
                        lists all the useful details similar to listobs in CASA or listr in AIPS.
  -t, --identify-targets
                        find targets for phasecal, science and bright cal for FF
  -r, --find-refant     find refant by checking TSYS info
  -s SPLIT_SOURCE, --split-source SPLIT_SOURCE
                        comma separated values of SOURCE_ID is used to select Sources to create a new FITS file.
  --to-B TO_B           Convert to B1920

Calibrator List:
  
                               use source based search on the calibrator list file.

  -S SEARCH_SOURCE, --search-source SEARCH_SOURCE
                        Comma separated list of sources to look for all the values
  -C COLUMN, --column COLUMN
                        Comma separated values of column to get value for
```

# Installation
After you clone the repo
```bash

$ cd vasco
$ pip install .
```


# Development
After you clone the repo
```bash

$ cd vasco
$ pip install -e .[dev]
```