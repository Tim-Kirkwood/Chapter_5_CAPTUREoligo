# Chapter_3_CAPTUREoligo
## Data preparation
Download your genome of interest and place in folder A.  If the genome consists of multiple contigs, place seperate genbank files corresponding to each contig into folder A.  

## Run code
Save CAPTUREoligo.py to folder A.  Assign the full path for folder A to the variable ```dir_with_input_gbks``` (line 557). Assign your preferred results folder name to the variable ```results_dir``` (line 553). ```concatenated_gbk_filename``` should be your preferred filname for the concatenated genbank file generated from the seperate genome genbank files in folder A.  Update the variables between lines 563 and 605 with preferred search spans, thermodynamic/GC/sequence limits, PAM characteristics etc - see code comments for full description.

Once CAPTUREoligo.py has finished running, you will have the following directory structure:
```
folder A/results_dir/forward_oligos
                    /reverse_oligos
```

The forward and reverse oligo folders will contain N subfolders, where N is the number of CAPTURE regions being designed.  Each subfolder will have the format first_cds_to_last_cds, where first_cds and last_cds are the locus tags defined in ```codon_boundary_list``` (line 595).  Each subfolder will contain the protospacer/gRNA/primer sequences and accompanying data in CSV and txt format, as well as a series of dynamic html plots showing the distributions of each possible filter parameter for the identified sequences.  This can be used to fine tune limits for the various design thresholds offered in the pipeline (secondary structure formation limits, etc). Data from the CSV files can then be fed into an oligo design software like snapgene if users wish to further prune the sequences (I used SnapGene to remove oligos that were less specific).   
