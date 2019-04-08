## species-energy-simulation


Code for conducting eco-evolutionary simulations of diversification and dispersal of species with 
and without energetic constraints as described in Hurlbert & Stegen, 2014a, "When should species 
richness be energy limited, and how would we know?" *Ecology Letters*. [DOI: 10.1111/ele.12240](http://onlinelibrary.wiley.com/doi/10.1111/ele.12240/abstract) and in Hurlbert & 
Stegen, 2014b, "On the processes generating latitudinal richness gradients: identifying diagnostic
patterns and predictions", *Frontiers in Genetics*. [http://dx.doi.org/10.3389/fgene.2014.00420](http://journal.frontiersin.org/Journal/10.3389/fgene.2014.00420/abstract)

## Setup
Requirements: R 3.0 or greater with the following packages installed and the following scripts sourced:

```R
package.vector = c('ape',
		   'caper',
		   'paleotree',
		   'plyr',
		   'foreach',
		   'doParallel')

install.packages(package.vector)

source('code/run_sim.r')
source('code/analyze_sim.r')
source('code/senc_sim_fun.r')
source('code/senc_analysis_fun.r')
source('code/supplemental_analysis_functions.r')
```

Cloning the Github repository will set up the following subdirectories:  

* code: most of the R scripts and functions  
* archived_sim_output: zipped simulation output for a subset of the sims conducted in the paper  
* raw_sim_output: folder in which new simulation output will be stored  
* analysis_output: folder in which new analysis output will be stored  


## Run simulations
The parameters governing a given simulation are specified within the SENC_Master_Simulation_Matrix.csv
file. These parameters include:

<table>
  <tr>
    <td>sim.id</td>
    <td>a unique numeric identifier referring to a set of simulation parameters</td>
  </tr>
  <tr>
    <td>status</td>
    <td>"completed" if the simulation has finished, or "to.run" if it is pending</td>
  </tr>
  <tr>
    <td>reg.of.origin</td>
    <td>region of origin, either "tropical" or "temperate"</td>
  </tr>
  <tr>
    <td>w</td>
    <td>governs the strength of environmental filtering, small values = stronger</td>
  </tr>
  <tr>
    <td>alpha</td>
    <td>per individual speciation probability</td>
  </tr>
  <tr>
    <td>beta</td>
    <td>per individual dispersal probability</td>
  </tr>
  <tr>
    <td>sigma_E</td>
    <td>governs the strength of niche conservatism, small values = stronger</td>
  </tr>
  <tr>
    <td>carry.cap</td>
    <td>"on" if limits exist to community abundance; else "off"</td>
  </tr>
  <tr>
    <td>energy.gradient</td>
    <td>"on" if carry.cap varies linearly across the gradient; else "off"</td>
  </tr>
  <tr>
    <td>max.K</td>
    <td>num of individuals that can be supported in region with the highest carrying capacity</td>
  </tr>
  <tr>
    <td>num.of.bins</td>
    <td>number of spatial bins</td>
  </tr>
  <tr>
    <td>max.time</td>
    <td>maximum number of time steps to run simulation</td>
  </tr>
  <tr>
    <td>max.richness</td>
    <td>maximum number of species that can accumulate before simulation breaks off</td>
  </tr>
  <tr>
    <td>replicate</td>
    <td>replicate ID if simulation is a repeat of an existing parameter combination</td>
  </tr>
  <tr>
    <td>disturb_frequency</td>
    <td>frequency of disturbance in number of time steps</td>
  </tr>
  <tr>
    <td>temperate_disturb_intensity</td>
    <td>fraction of individuals killed in disturbance event at temperate end of gradient</td>
  </tr>
  <tr>
    <td>tropical_disturb_intensity</td>
    <td>fraction of individuals killed in disturbance event at tropical end of gradient</td>
  </tr>
  <tr>
    <td>specn.gradient</td>
    <td>"on" if speciation rate (per individual probability) varies linearly across the gradient; else "off"</td>
  </tr>
  <tr>
    <td>specn.factor</td>
    <td>factor by which speciation rate varies from the highest to the lowest end of the gradient</td>
  </tr>
  <tr>
    <td>gamma</td>
    <td>The exponential decay rate describing how extinction probability decreases with population size</td>
  </tr>
</table>

The simulations reported on in the [*Ecology Letters*](http://onlinelibrary.wiley.com/doi/10.1111/ele.12240/abstract) paper correspond to the following sim.id's:

Zero sum energy gradient, temperate origin: c(4075:4084, 4275:4364)  
Zero sum energy gradient, tropical origin: c(4065:4074, 4185:4274)  
No zero sum constraint, temperate origin: c(3565:3664)  
No zero sum constraint, tropical origin: c(3465:3564)  

The simulations reported on in the [*Frontiers in Genetics*](http://journal.frontiersin.org/article/10.3389/fgene.2014.00420/abstract) paper correspond to the following sim.id's:

Energy gradient, temperate origin: c(4075:4084)  
Energy gradient, tropical origin: c(4065:4074)  
Pure niche conservatism, temperate origin: c(3565:3574)  
Pure niche conservatism, tropical origin: c(3465:3474)  
Speciation gradient, temperate origin: c(5535:5544)  
Speciation gradient, tropical origin: c(5525:5534)  
Disturbance gradient, temperate origin: c(5635:5644)  
Disturbance gradient, tropical origin: c(5625:5634)

(Note that phylogenies in Figure 2 were estimated at 30,000 timesteps rather than 100,000 for the Energy, Disturbance, and 
Speciation Gradient scenarios. For the Disturbance Gradient, this is analogous to examining sim 3865 instead of 5625. For Energy
and Speciation Gradient scenarios, a phylogeny of species extant at t=30,000 was created using code like this:
```
sim = output.unzip('raw_sim_output', simID)
all.pops = sim$all.populations
extant.pops = subset(all.pops, time.of.origin <= 30000 & time.of.extinction > 30000)
phy = sim$phylo.out
extant_phy = drop.tip(phy, phy$tip.label[!phy$tip.label %in% extant.pops$spp.name])
```
The extant phylogenies were then analyzed using [BAMM](http://bamm-project.org/) and [BAMMtools](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12199/abstract). See this [Github repo]
(https://github.com/ahhurlbert/bamm-simulations) for BAMM-related analysis code.)

-----

To run simulations with novel parameter combinations, add a line or lines specifying those parameter combos
to the SENC_Master_Simulation_Matrix.csv file.  

Then, as an example, if you'd like to run sim IDs 3325 to 3334 on your local machine, parallelizing over 2 processors, type
`run_sim(3325:3334, local = T, num_cores = 2)`

All simulation output will appear in the `raw_sim_output` folder. You can check on the status of a large simulation run by 
looking in the 'log.txt' file created in this folder.

The inputs required for running simulations and the output
files produced are listed below.

<table>
  <tr>
    <td>SIMULATION OUTPUT:</td>
    <td>(all saved in 'raw_sim_output' folder; XXX specifies simID)</td>
  </tr>
  <tr>
    <td>SENC_all.pops_simXXX.csv</td>
    <td>Simulation output containing all spatial, temporal and trait information for every  
    population to arise over the course of the simulation.</td>
  </tr>
  <tr>
    <td></td>
    <td>region:   spatial region number</td>
  </tr>
  <tr>
    <td></td>
    <td>spp.name:   species identifier</td>
  </tr>
  <tr>
    <td></td>
    <td>extant:   species extant as of the most recent time step (1) or extinct (0)</td>
  </tr>
  <tr>
    <td></td>
    <td>env.opt:   environmental optimum of the species, ranging from ~0-40</td>  
  </tr>
  <tr>
    <td></td>
    <td>time.of.origin:   time at which the species appeared in the region. If this is a first  
      appearance due to having newly speciated, then time.of.origin should equal  
      time.of.sp.origin. More likely, the appearance is due to colonization from an  
      adjacent region.</td>  
  </tr>
  <tr>
    <td></td>
    <td>time.of.extinction:   time at which the species disappears from the region. Most often  
      this will reflect the timing of a local extinction.</td>  
  </tr>
  <tr>
    <td></td>
    <td>time.of.sp.origin:   time at which the species first appeared in the simulation due to  
      speciation.</td>  
  </tr>
  <tr>
    <td></td>
    <td>time.of.sp.extinction:   time at which the species went globally extinct. If it did not  
      go extinct, then its value is the number of time steps the simulation was run plus 1.</td>  
  </tr>
  <tr>
    <td></td>
    <td>reg.env:   the average temperature of the region.</td>  
  </tr>
  <tr>
    <td></td>
    <td>carry.cap:   the carrying capacity (over all species) of the region.</td>
  </tr>
  <tr>
    <td>SENC_time.rich_simXXX.csv</td>
    <td>Summary of the number of species in each region at each time step in the simulation.</td>
  </tr>
  <tr>
    <td>SENC_phylo_simXXX.tre</td>
    <td>Object of class 'phylo' describing phylogenetic relationships among simulation species.</td>
  </tr>
  <tr>
    <td>SENC_params.out_simXXX.csv</td>
    <td>Simulation parameters used for simulation XXX.</td>
  </tr>
</table>


## Analyze simulation output
To analyze raw simulation output created using 'run_sim', use the 'analyze_sim' function. With this
function, you can specify whether you want to analyze patterns for the root clade only (i.e., for the
entire clade of species that diversified over the course of the simulation; root.only = 1) or for all
possible subclades nested within and including the root clade (root.only = 0; this will obviously
take much longer). You may also specify whether you want to analyze patterns at multiple points in 
time ("timeslices"), or just the pattern present at simulation's end (the default).

You can also specify whether to analyze some of the example simulation output that has already been
run (sim_dir = "archived") or simulation output that you have newly created (sim_dir = "new"). The
'archived_sim_output' folder contains simulation output for 10 of the 100 simulations for each of the
4 combinations of region of origin and zero sum condition. The entire simulation output (~7 GB) can be
found in the [Dryad Data Repository associated with this paper](http://datadryad.org/resource/doi:10.5061/dryad.4r8p3).

For example, to analyze simulation output for sim IDs 4065:4074 for all possible subclades but for
just the final timeslice on your local machine, parallelizing over 2 processors, type 
`analyze_sim(4065:4074, local = T, num_cores = 2, root.only = 0, num.of.time.slices = 1, sim_dir = "archived")`

The table below provides more details on the arguments to this function, and describes the
analysis output

<table>
    <td>ARGUMENTS</td><td></td>
  </tr>
  <tr>
    <td>which.sims</td>
    <td>simulation IDs to analyze</td>
  </tr>
  <tr>
    <td>local</td>
    <td>toggle for analyses on local machine (1) or computing cluster (0)</td>
  </tr>
  <tr>
    <td>num.cores</td>
    <td>number of processors available for parallel processing on a local machine</td>
  </tr>
  <tr>
    <td>root.only</td>
    <td>analyze patterns for the root clade only (1) or for all possible subclades including the root (0)</td>
  </tr>
  <tr>
    <td>num.of.time.slices</td>
    <td>number of points in time to conduct analyses; if 1, then analysis will be conducted at the final time step in the simulation</td>
  </tr>
  <tr>
    <td>which.time.slices</td>
    <td>a vector of points in time to analyze if those time points are irregularly spaced</td>
  </tr>
  <tr>
    <td>time.sequence</td>
    <td>a vector of points in time to analyze if those time points are regularly spaced; NB: JCS is responsible for these crazy time vectors</td>
  </tr>
  <tr>
    <td>min.num.spp</td>
    <td>the minimum number of species in a clade needed to proceed with analysis (i.e., clades smaller than this will be skipped)</td>
  </tr>
  <tr>
    <td>sim_dir</td>
    <td>repository directory where the simulation output you wish to analyze resides; "archived" if analyzing provided example simulation output, "new" if analyzing newly created output</td>
  </tr>
  <tr>
    <td></td><td></td>
  </tr>
  <tr>
    <td>ANALYSIS OUTPUT:</td><td></td>
  </tr>
  <tr>
    <td>SENC_Stats_simXXX.csv</td>
    <td>For each subclade and at
    each specified timeslice within simulation XXX, the following information is stored:
    simID, cladeID, time, correlation coefficients and p-values for correlations among
    richness, environment, region, MRD, and PSV, and Pybus & Harvey's gamma statistic.</td>
  </tr>
  <tr>
    <td>summary_output_simXXX.csv</td>
    <td>One row for each simulation with summary information regarding the simulation
    parameters, the total number of regions over which species occurred, the number of 
    extant species, the number of species going extinct over the simulation, the number
    of clades that were skipped in analysis because they were smaller than the min.num.spp
    threshold, the number of timeslices that were skipped because there were insufficient
    species (usually timeslices very early on under a limited range of parameters), and
    the total number of rows of output in that simulation's SENC_Stats_simXXX.csv file.</td>
  </tr>
</table>

## Notes on running and analyzing simulations on a High Performance Cluster
Clusters differ in their setup, and so adjustments to the code may be required.
As written, the code assumes that the cluster is configured for Message Passing Interface (MPI),
which handles the parallelization of the computing.

To run simulations on the UNC KillDevil cluster (which uses LSF to execute batch jobs), 
first clone the repository, and then from the repository home directory, type  
`$ bsub -o out.%J -n 100 -a openmpi mpirun Rscript run_sim_on_cluster.r 3765 3864`  

In this example, this command requests 100 nodes to run simulations 3765 to 3864. A job output file
is created called 'out.%J' where %J is the job number. Similarly, to analyze simulations on the cluster:

`$ bsub -o out.%J -n 100 -a openmpi mpirun Rscript run_analysis_on_cluster.r 3765 3864 30`  

In this case, an additional parameter at the end (here, 30) specifies the number of points in time over
which to analyze simulation output. Whether running simulations or conducting analyses, results will be 
stored in the `raw_sim_output` or `analysis_output` subdirectories.

Note: The analysis of large phylogenetic trees is memory-intensive. You may need to increase the memory allocation on the cluster using the -M option, in which you can set the memory limit in GB (default is 4):

`$ bsub -M 8 -o out.%J -n 100 -a openmpi mpirun Rscript run_analysis_on_cluster.r 3765 3864 30`


## Duplicating manuscript figures from Hurlbert & Stegen 2014a, *Ecology Letters*:

### Figures 2 and 4
metrics.through.time.r

### Figure 3
KT_scenarios_corr_plots_plus_sebastes.r

### Figure S1
lat.grad.time.plot.r

### Figure S2
metrics.through.time.sigma.variation.r

### Figure S3
pacific_npp.r  
sebastes.r
