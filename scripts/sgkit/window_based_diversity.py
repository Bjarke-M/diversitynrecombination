import sgkit as sg
import numpy as np
from sgkit.io.vcf import vcf_to_zarr
#import bcfs
#vcf_to_zarr("/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/liftover/Daubentonia_madagascariensis/bcfs/Daubentonia_madagascariensis_lifted.bcf.gz", "output.zarr")
#vcf_to_zarr("/home/bjarkemp/primatediversity/people/bjarkemp/diversitynrecombination/data/liftover/Gorilla_gorilla_gorilla/bcfs/Gorilla_gorilla_gorilla_lifted.bcf.gz", "output_gorilla.zarr")
ds = sg.load_dataset("output_gorilla.zarr")
# make windows in the data frame:
print(ds)
#test = sg.individual_heterozygosity(ds)['call_heterozygosity']
#print(test)