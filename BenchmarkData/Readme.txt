BenchmarkData contains the simulated data analyzed in the Rensing lab without any knowledge about the underlying truth. The configuration parameters had to be chosen as good as possible like in real application settings without this knolwedge

GroundTruth contains the data AND information about differential regulation:
methylated		1,2,3,4	  methylation state (termed S_1, S_2, S_3, S_4 in our publication [2])
island_state		1,2	CpG island state (2) or deserts (1)
diff_methylated		Magnitude of the methylation difference (termed as $\Delta$ in our publication [2]).

The meaning of these parameter are originate from the WGBSS_Suite [1] used for simulation.


For unbiased future benchmark studies, please strictly prohibit knowlegde about the truth.

References:

[1] Rackham, O. J. L., Dellaportas, P., Petretto, E., and Bottolo, L. (2015). WGBSSuite: simulating whole-genome bisulphite sequencing data and benchmarking
differential DNA methylation analysis tools. Bioinformatics, 31(14), 2371–2373

[2] Kreutz et al., submitted.

