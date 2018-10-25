
# PGSMs
Rprof("profile1.out", line.profiling=TRUE)
eval(parse(file = "C:/Users/dt915355/Dropbox/Postdoc/Projects/ParticleGibbsforSBMs/Code/Dashboard - PGSMs.R", 
           keep.source=TRUE))
Rprof(NULL)
summaryRprof("profile1.out", lines = "show")$by.self


# Gibbs
Rprof("profile1.out", line.profiling=TRUE)
eval(parse(file = "C:/Users/dt915355/Dropbox/Postdoc/Projects/ParticleGibbsforSBMs/Code/Dashboard - Gibbs.R", 
           keep.source=TRUE))
Rprof(NULL)
summaryRprof("profile1.out", lines = "show")$by.self

# PGSMs & Gibbs
Rprof("profile1.out", line.profiling=TRUE)
eval(parse(file = "C:/Users/dt915355/Dropbox/Postdoc/Projects/ParticleGibbsforSBMs/Code/Dashboard - Main.R", 
           keep.source=TRUE))
Rprof(NULL)
summaryRprof("profile1.out", lines = "show")$by.self


