# UMI_vis
UMI visualization for Fulcrum Genomics UMI, who logs the UMI sequence in RX tag in BAM format.

## Install
```
git clone git@github.com:dazhouze/UMI_vis.git
```

## Run
umi_vis.py have 2 command, dis and vis. "dis" is for UMI distribution analysis and "vis" is for UMI visualzation in given. coordinate.
```
# for details please run
python3 umi_vis.py 
```
The "vis" command will generate SVG plot for given genome region and reads with one UMI will be colored grey and more than 1 UMI will be colored randomly except grey.
