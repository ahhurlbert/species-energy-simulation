##species-energy-simulation


Code for conducting eco-evolutionary simulations of diversification and dispersal of species with and without energetic constraints as described in Hurlbert &
Stegen, *Ecology Letters*.

#Setup
Requirements: R 3.0 with the following packages installed:

```sh
package.vector = c('ape','permute','nlme','vegan','picante','mvtnorm','caper','paleotree','plyr','phytools','apTreeshape','foreach','doParallel')

install.packages(package.vector)
```

#Run simulations
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


The simulations reported on in the paper correspond to the following sim.id's:

Zero sum energy gradient, temperate origin: c(4075:4084, 4275:4364)  
Zero sum energy gradient, tropical origin: c(4065:4074, 4185:4274)  
No zero sum constraint, temperate origin: c(3565:3664)
No zero sum constraint, tropical origin: c(3465:3564)

To run simulations with novel parameter combinations, add a line or lines specifying those parameter combos
to the SENC_Master_Simulation_Matrix.csv file. 

FILE:   sims_only_workflow_local.r
FILE INPUTS:   
 DATA
  SENC_Master_Simulation_Matrix.csv 
    Table listing the parameter combinations of every simulation to be run.
 SCRIPTS
  senc_sim_fun.r
    Runs an entire simulation for one sim ID (i.e., one combination of simulation parameters)

VARIABLE INPUTS:
  which.sims
    The sim IDs for which simulations should be run

FILE OUTPUTS:   (all saved in 'raw_sim_output' folder; XXX specifies simID)
  SENC_all.pops_simXXX.csv
    Simulation output containing all spatial, temporal and trait information for every
    population to arise over the course of the simulation.
    region:   spatial region number
    spp.name:   species identifier
    extant:   species extant as of the most recent time step (1) or extinct (0)
    env.opt:   environmental optimum of the species, ranging from ~0-40
    time.of.origin:   time at which the species appeared in the region. If this is a first
      appearance due to having newly speciated, then time.of.origin should equal
      time.of.sp.origin. More likely, the appearance is due to colonization from an
      adjacent region.
    time.of.extinction:   time at which the species disappears from the region. Most often
      this will reflect the timing of a local extinction.
    time.of.sp.origin:   time at which the species first appeared in the simulation due to
      speciation.
    time.of.sp.extinction:   time at which the species went globally extinct. If it did not
      go extinct, then its value is the number of time steps the simulation was run plus 1.
    reg.env:   the average temperature of the region.
    carry.cap:   the carrying capacity (over all species) of the region.
  SENC_time.rich_simXXX.csv
    Summary of the number of species in each region at each time step in the simulation.
  SENC_phylo_simXXX.tre
    Object of class 'phylo' describing phylogenetic relationships among simulation species.
  SENC_params.out_simXXX.csv
    Simulation parameters used for simulation XXX.


#Analyze simulation output
The following code creates derived data and statistics from the raw simulation output.

FILE:   analyzing_sims.r
---------------------------------
FILE INPUTS:   
 DATA
  SENC_Master_Simulation_Matrix.csv 
    Table listing the parameter combinations of every simulation to be run.
  senc.out.XXX.zip
    Zipped output files where XXX is the simulation id, ranging from 1 to ~4000
  (OR UNZIPPED:)
  SENC_all.pops_simXXX.csv
  SENC_time.rich_simXXX.csv
  SENC_phylo_simXXX.tre
  SENC_params.out_simXXX.csv
 SCRIPTS
  analysis_workflow_Olympus.r
    complete analysis of a single simulation
  unzipping_files.r
    Unzips simulation output of the form senc.out.XXX.zip, stores info in sim.results
  reg_calc_and_analysis.r
    
  extinct.calc.r
    Calculates an extinction rate for each region
  make.phylo.jimmy.fun.r
    Makes a phylogeny from the edge and edge.length info collected during the simulation
  lat.grad.time.plot.r
    Makes the lat_grad_thru_time_simXXX.pdf plots.
  clade.origin.corr.plot.r
    Makes the corr_vs_cladeage_simXXX.pdf plots.
  clade.exmpl.figs.r 
    Makes clade_example_figs_simXXX.pdf plots. Currently broken.

VARIABLE INPUTS:   
  which.sims
    simulation IDs to analyze
  local
    toggle for analyses on local machine (1) or computing cluster (0)
  num.cores
    number of processors available for parallel processing on a local machine


FILE OUTPUTS:
  SENC_Stats_simXXX.csv
    Often referred to as the "stats" file or "output" file. For each subclade and at
    each specified timeslice within simulation XXX, the following information is stored:
    simID, cladeID, time, correlation coefficients and p-values for correlations among
    richness, environment, region, MRD, and PSV, Pybus & Harvey's gamma statistic, and
    Blomberg's K for both environmental optimum and mean region of occurrence.
  summary_output_simXXX.csv
    One row for each simulation with summary information regarding the simulation
    parameters, the total number of regions over which species occurred, the number of 
    extant species, the number of species going extinct over the simulation, the number
    of clades that were skipped in analysis because they were smaller than the min.num.spp
    threshold, the number of timeslices that were skipped because there were insufficient
    species (usually timeslices very early on under a limited range of parameters), and
    the total number of rows of output in that simulation's SENC_Stats_simXXX.csv file.