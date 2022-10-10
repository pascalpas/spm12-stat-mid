# spm12-stat-mid

Preprocessing, QC and statistical analysis of the Monetary Incentive Delay task (MID) 

SPM12 required

## Usage
```
install and add spm12 to MATLAB path
add spm12-stat-mid and subfolders to path, i.e. addpath(genpath('spm-stat-mid')
```

### Preprocessing

> preproc_2DEPI(subjectDir,subjectName,anatomy,func,TR,mbFactor)
```
Input:
subjectDir:     main subjectDir
subjectName:    subject name (beginning of filenames)
anatomy:        anatomical dir name
func:           functional dir name
TR:             TR in seconds
mbFactor:       multiband factor for slice timing  (e.g. 3)
```

### Statistical analysis

> SPM_stat_job_MID(subjectDir,funcDir,behavDir,TR,outputDir)
```
Input:
subjectDir:     main subjectDir
funcDir:        directory containing preprocessed functionals
behavDir:       dir containing behavfile
TR:             TR in seconds
OutputDir:      Output directory for statistical results

```
