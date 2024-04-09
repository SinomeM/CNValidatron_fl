# CNValidatron/UKB CNVs TO-DO list

## Training and validation

- 3 classes instead of 5
- ADD negative examples. Regions with no CNV at all and locally noisy samples
- minimu number of SNPs in the image? (avoid too sparse pictures).
  carefull if one side is empty
- on the oppostite side, downsample SNPs when the CNV is very large.
  The general idea is to make small and large CNV more similar
- In general, it might be nice to have the pixel density more or less
  constant across the image
- Now that the image is larger the simple model might be too small,
  test whether making it larger helps (expecially on the combined training set, UKB + deCODE)


## Binning

-lorem


## Association analysis
- make the center position of each bin the "POS" columns and treat them as
  clumped SNPs???
- Especially when using large bin, it could be considered computing the
  IOU (super easy from the long table). Maybe bin size should be set at
  a point were the average IOU is quite high (meaning CNVs overlapping
  a bin tend to overlap it fully)
- Regarding the previous, Morten suggested there are plink format
  that support fractional values, would be nice to use those and
  for example compute the actual proportion of the bin covered by
  any CNVs, this way I can use larger bins without sacrificing too
  much specificity


## CNVRs

- Test whether we can use the IOU directly as connection, rather than
  0/1 for those connection above a certain IOU. I suspect it will
  require more RAM of that available on my laptop
- Doing the above but still using a IOU filter on which CNV to consider
  should not have any impact on memory use


## Misc

- The stichting function must be a lot faster
- The real number of SNPs (after stitching I mean) can be computed with
  the PNGs and used after the model has run (for filtering)
- There is a hmm specifically for affy in PennCNV, might be worth re
  running the calling
- not doing so might raise some questions (especially on the low number
  of CNV call per sample)
- for the training I don't think it's a problem, might be for the paper
- Need to ask Jesper to do it
- deCODE trios are very good for positive examples (True dels and dups),
  ~95% predicted true by genealogy above 30 SNPs are true
- For the paper, compare the efficacy of CNValidatron vs different set
  of cutoffs from "new quality measure" etc paper


## Notes

### UKB calibration
- if a sample have multiple CNVs are they more likely to occur close by
  than random?
- duplicated lines in the calibration set ??? (too few cnvs in the
  lower bins??) (edited)
- there are at least two large "predicted true deletions" that are
  clearly false in the validation set
- pred_prob is not properly calibrated but it's fine. If one set the
  minimum at 0.9 can expect "at lesst 90% of the predicted true are
  actually true, most likely way more"
