# BELCovAge: COVID-19 Age-Related Heterogeneity Analysis in Belgium

This repository hosts the codebase for replicating the findings of a longitudinal study examining age-related heterogeneity in COVID-19 transmission across Belgium. Detailed in the preprint manuscript below, this study highlights the role of different age groups during the pandemic.

**Citation:**
Angeli, L., Caetano, C., Franco, N., Coletti, P., Faes, C., Molenberghs, G., Beutels, P., Abrams, S., Willem, L., and Hens, N., 2024. _Insights into the role of children in the COVID-19 pandemic in Belgium: a longitudinal sensitivity analysis._ [DOI: 10.21203/rs.3.rs-4324206/v1](https://doi.org/10.21203/rs.3.rs-4324206/v1)


## Overview
This repository contains the R scripts and data necessary to replicate the findings of our longitudinal study on age-related heterogeneity in COVID-19 transmission across Belgium. This research is detailed in the preprint manuscript listed below, emphasizing different age groups' distinct roles during the pandemic.

## System & Software Requirements
- **R Version**: 4.3.0
- **RStudio Version**: 2023.12.0.369

These versions were used for the analysis, ensuring the reproducibility of the results. Please install these versions or later to avoid compatibility issues.

## Repository Contents

### Folders:
- `data`: Contains the data needed to reproduce the analysis.
- `fatigue_correction`: Includes estimated effects of fatigue when filling out CoMix social contact surveys [3.].
- `R`, `R_socialmix`: R scripts adapted from the SOCRATES project[1.] to handle Comix social contact data[2.].

### Scripts:
- `add_functions.R`: Defines necessary functions for the analysis.
- `helpers.R`: Contains the necessary functions to run supplementary analysis for the study.
- `preliminary.R`: Installs necessary packages and runs preliminary operations.
- `wave_specific_an.R`: Processes waves of the CoMix survey and generates variables required for the analysis.
- `long_analysis_BE.R`: Generates the results of the analysis.

## Citations
1. Willem, L., Van Hoang, T., Funk, S., et al. SOCRATES: an online tool leveraging a social contact data-sharing initiative to assess mitigation strategies for COVID-19. BMC Res Notes 13, 293 (2020). [https://doi.org/10.1186/s13104-020-05136-9](https://doi.org/10.1186/s13104-020-05136-9)
2. CoMix social contact data (Belgium), Coletti, Pietro; Wambua, James; Gimma, Amy; Willem, Lander; Vercruysse, Sarah; Vanhoutte, Bieke; Van Zandvoort, Kevin; Edmunds, John; Beutels, Philippe; Hens, Niel [https://doi.org/10.5281/zenodo.4035001](https://doi.org/10.5281/zenodo.4035001)
3. Loedy, N., Coletti, P., Wambua, J., et al. Longitudinal social contact data analysis: insights from 2 years of data collection in Belgium during the COVID-19 pandemic. BMC Public Health 23, 1298 (2023). [https://doi.org/10.1186/s12889-023-16193-7](https://doi.org/10.1186/s12889-023-16193-7)
4. Angeli, L., Caetano, C., Franco, N., Coletti, P., Faes, C., Molenberghs, G., Beutels, P., Abrams, S., Willem, L., and Hens, N., 2024. Insights into the role of children in the COVID-19 pandemic in Belgium: a longitudinal sensitivity analysis. [https://doi.org/10.21203/rs.3.rs-4324206/v1](https://doi.org/10.21203/rs.3.rs-4324206/v1)
5. Abrams, S., Wambua, J., Santermans, E., Willem, L., Kuylen, E., Coletti, P., Libin, P., Faes, C., Petrof, O., Herzog, S.A. and Beutels, P. Modelling the early phase of the Belgian COVID-19 epidemic using a stochastic compartmental model and studying its implied future trajectories. Epidemics, 35, p.100449(2021).
[https://doi.org/10.1016/j.epidem.2021.100449]
6. Willem, L., Abrams, S., Franco, N., Coletti, P., Libin, P.J., Wambua, J., Couvreur, S., André, E., Wenseleers, T., Mao, Z. and Torneri, A. The impact of quality-adjusted life years on evaluating COVID-19 mitigation strategies: lessons from age-specific vaccination roll-out and variants of concern in Belgium (2020-2022). BMC public health, 24(1), p.1171(2024).
[https://doi.org/10.1186/s12889-024-18576-w]

## How to Use
To replicate the study results:
1. Clone this repository.
2. Ensure R and RStudio are installed with the versions specified above.
3. Run `preliminary.R` to install required packages.
4. Execute `wave_specific_an.R` and `long_analysis_BE.R` sequentially to generate the results and plots as described in the manuscript.
5. If you are just interested in the final results, it is sufficient to run `long_analysis_BE.R`, which will automatically run all the dependencies.
6. Expected run time: 13 mins
## Contact
For issues, questions, or contributions, please open an issue on this GitHub repository or contact the repository maintainers directly.
