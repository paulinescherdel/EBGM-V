This code repository is associated with the publication "Development, internal and external evaluation of an artificial intelligence algorithm for child growth monitoring". 

We provide the codes : 
- to model individual height trajectories using the Jenss-Bayley non-linear mixed model according to sex, 
- to develop age-specific predictive models by age range and obtain thresholds for a pre-defined specificity.

Fictive samples of the dataset are included in this repository as sample_G.rda (for boys) and sample_F (for girls). 
The full datasets are not publicly released due to privacy and licensing constraints. 

These samples include as variables: 
- id
- sex (1=boy, 2=girl)
- age (in days)
- height (in cm)
- x (outcome, 1=disease, 0=non-disease)
- ts (Turner syndrome, 1=oui, 0=non)
- ghd (Growth Hormone Deficiency, 1=oui, 0=non)
- referents (Non-disease, 1=oui, 0=non)
- controls (Modeling sample, 1=oui, 0=non)
