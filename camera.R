library(CAMERA)

setwd('/Rdata')
xset = readRDS('xcms_centroid_test_brain.Rdata')

# Create annotated object
an_xset   =  xsAnnotate(xset)

# Group after RT value of the xcms grouped peak
an_group  =  groupFWHM(an_xset, perfwhm = 0.6)

# Annotate Isotopes
an_iso    = findIsotopes(an_group, mzabs = 0.01) 

# Verify grouping
an_corr    =  groupCorr(an_iso, cor_eic_th = 0.75)

# Annotate Adducts
an_add   =  findAdducts(an_corr, polarity="positive")


peaklist =  getPeaklist(an_add) 
write.csv(peaklist, file='xs_annotated_test_brain.csv')
