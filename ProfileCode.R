

Rprof("profile1.out", line.profiling=TRUE)
eval(parse(file = "C:/Users/dt915355/Dropbox/Postdoc/Code/ParticleGibbsforSBMs/PerformPGSMsForSBM.R", 
           keep.source=TRUE))
Rprof(NULL)
summaryRprof("profile1.out", lines = "show")$by.self

