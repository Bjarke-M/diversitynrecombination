Overview of the data that will be created duing the project.

The folder will be devided into the following structure:

Data/ ->    
reference genome specific folders/: containing individual folders of the individuals/species that where mapped to this reference genome. ->
    LiftOver -> contain the lift over VCF and bcf files 
    Mask -> contains the liftOver BED files per individual, these are used to make a callability mask (they are generated on the gvcf files)
    
    Species_specific_bcfs -> reference assembly to whcih every species is mapped -> species -> nonpar: all the chromosome that is not par ->    females: alle the females in the sample, 
                                                                                                                                                males_no_x: all males, X chromosome not included in the bcf and a txt with the id for the males
                                                                                                                                                merged_non_male_X: merged bcf with females and males where the x from males have been removed
                                                                                                                                                {species}_nonpar.bcfs: isolated bcfs without par but including males and females and the male X
                                                                                            ->  original:   isolated species bcfs including par and male x
                                                                                            -> par: bcfs with only the PAR region
                                                                                            -> {species}_isolation.txt: a txt file of the species ids to that specific species used to filter the reference assembly multispecies bcf




There is also a metadata file, obtained from Erik S. for which i'll be ever greatful.


