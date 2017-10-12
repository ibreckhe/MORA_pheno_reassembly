# MORA_pheno_reassembly
Code and data suppporting "Theobald, Breckheimer, and HilleRisLambers. 2017. Climate drives phenological reassembly of a mountain wildflower meadow community. In press. Ecology." 

## References
Paper and Supplemental: doi:10.1002/ecy.1996 ([link](http://onlinelibrary.wiley.com/doi/10.1002/ecy.1996/full))

Code and Data: [![DOI](https://zenodo.org/badge/105659951.svg)](https://zenodo.org/badge/latestdoi/105659951)

**Abstract**

Spatial community reassembly driven by changes in species abundances or habitat occupancy is a well-documented response to anthropogenic global change, but communities can also reassemble temporally if the environment drives differential shifts in the timing of life events across community members. Much like spatial community reassembly, temporal reassembly could be particularly important when critical species interactions are temporally concentrated (e.g. plant-pollinator dynamics during flowering). Previous studies have documented species-specific shifts in phenology driven by climate change, implying that temporal reassembly, a process we term “phenological reassembly,” is likely. However, few studies have documented changes in the temporal co-occurrence of community members driven by environmental change, likely because few datasets of entire communities exist. We addressed this gap by quantifying the relationship between flowering phenology and climate for 48 co-occurring subalpine wildflower species at Mount Rainier (Washington, USA) in a large network of plots distributed across Mt. Rainier’s steep environmental gradients; large spatio-temporal variability in climate over the six years of our study (including the earliest and latest snowmelt year on record) provided robust estimates of climate-phenology relationships for individual species. We used these relationships to examine changes to community co-flowering composition driven by ‘climate change analog’ conditions experienced at our sites in 2015. We found that both the timing and duration of flowering of focal species was strongly sensitive to multiple climatic factors (snowmelt, temperature, and soil moisture). Some consistent responses emerged, including earlier snowmelt and warmer growing seasons driving flowering phenology earlier for all focal species. However, variation among species in their phenological sensitivities to these climate drivers was large enough that phenological reassembly occurred in the climate change analog conditions of 2015. An unexpected driver of phenological reassembly was fine-scale variation in the direction and magnitude of climatic change, causing phenological reassembly to be most apparent early and late in the season and in topographic locations where snow duration was shortest (i.e., at low elevations and on ridges in the landscape). Because phenological reassembly may have implications for many types of ecological interactions, failing to monitor community-level repercussions of species-specific phenological shifts could underestimate climate change impacts.

Please see the paper and supplemental material for details about the study design and methodology.

## Analysis Requirements

Code to reproduce the data analysis can be found in the `./code` directory.

To reproduce this work, you will need to install a suitable version (3.1 or later) of the [R computing environment](https://cran.r-project.org/), as well as the program [JAGS](http://mcmc-jags.sourceforge.net/) for performing MCMC estimation.

In addition to base R, you will also need to install the following packages and their dependencies:
* `rjags`
* `ggmcmc`
* `ggplot2`
* `reshape2`
* `dplyr`
* `gridExtra`
* `stringr`
* `ROCR`

The analysis was performed on a computer running MacOS Sierra with 8GB RAM. The code has not been tested on other platforms.

## Data and Metadata

Original data and metadata contributing to the analysis can be found in the `./data` directory.

## License

This work is shared under the [MIT License](https://www.tldrlegal.com/l/mit).

## Contact the Author

For questions, or to report any problems, please contact the authors Ian Breckheimer (`ibreckhe@u.washington.edu`), Elli Theobald (`ellij@u.washington.edu`), and Janneke Hille Ris Lambers (`jhrl@u.washington.edu`).

