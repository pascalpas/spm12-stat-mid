# MID_SPM12_stat

Preprocessing, QC and statistical analysis of the Monetary Incentive Delay task (MID) 

SPM12 required

## Usage

### Preprocessing

> preproc_2DEPI(subjectDir,subjectName,anatomy,func,scanTR,mbFactor)
```
Input:
subjectDir:     directoryname of subject;
subjectName:    subject name;
anatomy:        anatomical dir name
func:           functional dir name
scanTR:         TR of functionals (e.g. 1)
mbFactor:       multiband factor for slice timing  (e.g. 3)
```

### Statistical analysis

> SPM_stat_job_MID(analysisDir,funcDir,behavDir,TR,outputDir)
```
Input:
analysisDir:    main subjectDir
funcDir:        directory containing preprocessed functionals
behavFile:      path to behavioral file
TR:             TR in seconds
OutputDir:      Output directory for statistical results

```
