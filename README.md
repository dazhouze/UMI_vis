# UMI_vis
UMI visualization for Fulcrum Genomics UMI, who logs the UMI sequence in RX tag in BAM format.

## Install
```
git clone git@github.com:dazhouze/UMI_vis.git
```

## Run
umi_vis.py have 2 command, dis and vis. "dis" is for UMI distribution analysis and "vis" is for UMI visualzation in given. coordinate.
```
python3 umi_vis.py <dis/vis> <chr:start-end> <output> <bam1> [bam2 ...]
```
