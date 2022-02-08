This folder is meant to contain all resources necessary for running the workflow, for example reference sequences or databases.
Wherever feasible, they can also be downloaded programmatically via rules defined in the pipeline.

## qual_profile.txt

By default, NGmerge uses hard-coded profiles when determining the quality scores of overlapping bases.  There are separate profiles for cases where the R1 base and the R2 base match, and for when they do not match.  Those who do not wish to use these profiles have two alternative options:

```
  -w  <file>       Use given error profile for merged qual scores
```

The here provided qual_profile.txt is an adjusted version of the original quality profile for Phred+33 scores of Illumina 1.8+, which ranges from 0 to 41 instead of 0 to 40 in earlier versions. The only modification of the files is that the last column was copied and appended as an additional column.
