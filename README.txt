Outcome-informed reweighting of the CDC Social Vulnerability Index improves census-tract diabetes prediction in California.

CDC's SVI weights its 16 social vulnerability variables equally. This is a deliberate choice for general-purpose use, but one that sacrifices accuracy for any specific health outcome. We reweighted the same 16 variables against CDC PLACES diabetes prevalence using LASSO, then tested on held-out California census tracts.

Result: R² 0.70 vs. 0.47 for standard SVI (ΔR² = +0.24, 95% CI 0.22–0.27). The advantage holds under 10-fold CV, spline/decile flexibility for SVI, small-tract exclusion, and in 32 of 34 counties.

Why it matters: Health departments using SVI to target diabetes resources may be mis-ranking tracts, especially at the high-risk tail where allocation decisions happen.

Data: CDC/ATSDR SVI 2022 (CA) + CDC PLACES 2023