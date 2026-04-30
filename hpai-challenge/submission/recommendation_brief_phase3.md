# Recommendation brief: HPAI epidemic — final assessment

## Epidemic outcome

The HPAI H5N0 epidemic on Jolly Island ran from 22 December 2025 to 7 April 2026, producing 560 confirmed outbreaks across 107 days. Chicken farms accounted for 75% of cases (broiler\_2 alone: 46%). Incidence peaked in late January at approximately 18 cases per day, then declined steadily through February and March. The last case was confirmed on 7 April. Preventive culling removed 691 farms from the susceptible population, of which 74 turned out to be already infected.

## Restocking

Our fitted model does not identify a safe restocking date within the assessment window. The estimated spillover decay rate implies wildlife-mediated introduction risk persisting for months after the last case. This is probably an artefact of the model's spillover profile (a smooth exponential decay) rather than the true situation: the migration period has ended, and no cases have been reported for over a week.

**Practical recommendation.** If active surveillance detects no new infections through May, restocking could begin in June, prioritising farms furthest from former outbreak clusters. If the wildlife source has genuinely ceased, the risk from farm-to-farm transmission alone is negligible — all known infection chains have been broken.

## Preventive culling

We estimate that preventive culling averted approximately 130 outbreaks (95% CrI: −195 to 503) over the course of the epidemic. The wide interval reflects the difficulty of separating the effect of culling from other factors (susceptible depletion, confinement, declining spillover). In some model trajectories, removing culling leads to faster local burnout, producing a similar total. On balance, preventive culling reduced the overall case count, though its marginal benefit diminished as the epidemic declined.

## Culling capacity

At no tested capacity level (5 to 100+ farms per day) could the model achieve 95% probability of eradication by 13 February 2026. The reason is that culling addresses farm-to-farm transmission but cannot prevent new introductions from the wildlife reservoir. As long as environmental spillover continues, new outbreaks will occur regardless of how quickly detected farms are depopulated. Eradication required the natural end of the wildlife migration period, which occurred in late February–March.

## Key uncertainties

The spillover decay rate (δ = 0.018 day⁻¹, half-life ~39 days) is the most consequential parameter for the restocking and culling capacity questions. The model cannot distinguish a slowly decaying spillover process from the abrupt end of migration that occurred in reality. The confinement effect (50% susceptibility reduction) is assumed, not estimated. The counterfactual estimates (Q3, Q4) depend on the generative model approximating the observed epidemic dynamics, which it does within broad uncertainty bounds.
