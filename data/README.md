Overview of the data that will be created duing the project.

The folder will be devided into the following structure:

Data/ ->    
reference genome specific folders/: containing individual folders of the individuals/species that where mapped to this reference genome. ->
    LiftOver -> contain the lift over VCF and bcf files 
    Mask -> contains the liftOver BED files per individual, these are used to make a callability mask (they are generated on the gvcf files)
    


There is also a metadata file, obtained from Erik S. for which i'll be ever greatful.
