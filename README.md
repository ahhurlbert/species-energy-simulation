##species-energy-simulation


Code for conducting eco-evolutionary simulations of diversification and dispersal of species with 
and without energetic constraints as described in Hurlbert & Stegen, *Ecology Letters*.

##Setup
Requirements: R 3.0 or greater with the following packages installed and the following scripts sourced:

```R
package.vector = c('mnormt',
		   'rgl',
		   'ape',
		   'permute',
		   'nlme',
		   'vegan',
		   'picante',
		   'mvtnorm',
		   'caper',
		   'paleotree',
		   'plyr',
		   'phytools',
		   'apTreeshape',
		   'foreach',
		   'doParallel')

install.packages(package.vector)

source('code/run_sim.r')
source('code/analyze_sim.r')
source('code/senc_sim_fun.r')
source('code/senc_analysis_fun.r')
source('code/supplemental_analysis_functions.r')
```

##Run simulations
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
    <td>min.K</td>
    <td>num of individuals that can be supported in region with the lowest carrying capacity</td>
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
    <td>temperate_disturb_intensity</td>
    <td>fraction of individuals killed in disturbance event at temperate end of gradient</td>
  </tr>
  <tr>
    <td>tropical_disturb_intensity</td>
    <td>fraction of individuals killed in disturbance event at tropical end of gradient</td>
  </tr>
  <tr>
    <td>disturb_frequency</td>
    <td>frequency of disturbance in number of time steps</td>
  </tr>
</table>

The simulations reported on in the paper correspond to the following sim.id's:

Zero sum energy gradient, temperate origin: c(4075:4084, 4275:4364)  
Zero sum energy gradient, tropical origin: c(4065:4074, 4185:4274)  
No zero sum constraint, temperate origin: c(3565:3664)  
No zero sum constraint, tropical origin: c(3465:3564)  

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
    <td colspan="2">FILE: run_sim.r</td>
  </tr>
  <tr>
    <td>FILE INPUTS:</td><td></td>
  </tr>
  <tr>
    <td>SENC_Master_Simulation_Matrix.csv</td>
    <td>Table listing the parameter combinations of every simulation to be run.</td>
  </tr>
  <tr>
    <td>senc_sim_fun.r</td>
    <td>Runs an entire simulation for one sim ID (i.e., one combination of simulation parameters)</td>
  </tr>
  <tr>
    <td></td><td></td>
  </tr>
  <tr>
    <td></td><td></td>
  </tr>
  <tr>
    <td>FILE OUTPUTS:</td>
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


##Analyze simulation output
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
found in the Dryad Data Repository associated with this paper [insert link here].

For example, to analyze simulation output for sim IDs 4065:4074 for all possible subclades but for
just the final timeslice on your local machine, parallelizing over 2 processors, type 
`analyze_sim(4065:4074, local = T, num_cores = 2, root.only = 0, num.of.time.slices = 1, sim_dir = "archived")`

The table below provides more details on the arguments to this function, and describes the
analysis output

<table>
  <tr>
    <td colspan="2">FILE:   analyze_sim.r</td>
  </tr>
  <tr>
    <td>DATA INPUT REQUIRED</td><td></td>
  </tr>
  <tr>
    <td>senc.out.XXX.zip</td>
    <td>Zipped output files where XXX is the simulation id, ranging from 1 to ~4000</td>
  </tr>
  <tr>
    <td>(OR UNZIPPED:)</td><td></td>
  </tr>
  <tr>
    <td>SENC_all.pops_simXXX.csv</td>
    <td>raw data of species distributions, time of origin and extinction, and traits over the entire simulation</td>
  </tr>
  <tr>
    <td>SENC_time.rich_simXXX.csv</td>
    <td>summary of the number of species per region in each time step of the simulation</td>
  </tr>
  <tr>
    <td>SENC_phylo_simXXX.tre</td>
    <td>a phylogeny of all species (extinct and extant) generated by the simulation; can be read in R using read.tree() in 'ape'</td>
  </tr>
  <tr>
    <td>SENC_params.out_simXXX.csv</td>
    <td>the parameters associated with that particular simulation run</td>
  </tr>
  <tr>
    <td></td><td></td>
  </tr>
  <tr>
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
    <td>FILE OUTPUTS:</td><td></td>
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

##Duplicating manuscript figures

###Figures 2 and 4
metrics.through.time.r

###Figure 3
KT_scenarios_corr_plots_plus_sebastes.r

###Figure S1
lat.grad.time.plot.r

###Figure S2
metrics.through.time.sigma.variation.r

###Figure S3
pacific_npp.r  
sebastes.r
