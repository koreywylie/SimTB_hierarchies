3/17/2023 --kw--
  Customized scripts for creating simulated fMRI datasets with spatial hierarchies, for use with simtb (http://mialab.mrn.org/software).

Instructions:
  Copy below scripts to experiment directory.
  Note that sripts below must be in matlab path, with higher priority than simtb folder.


Novel scripts:
  drchrnd.m     : samples from Dirichlet distribution
  get_fractal.m : creates binary spatial maps with hierarchically nested fractal distribution

Modified Simtb scripts:
  simtb_create_sP.m : includes input for hierarchy parameters 
  simtb_makeSM.m    : includes input for hierarchy parameters
  simtb_SMsource.m  : includes generic spatial maps needed to construct hierarchy

Experiments:
  simtb_exp_hierarchy_5Levels_2Degree.m
  simtb_exp_hierarchy_3Levels_3Degree.m
  simtb_exp_hierarchy_3Levels_3Degree_Dirichlet.m