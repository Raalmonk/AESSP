# Limitations

AESSP is a screening-prioritization system, not a biological truth engine.

## Evidence limitations

- Abstracts and titles often omit culture conditions, assay details, and product-characterization methods.
- Reported yield, molecular weight, and branching values may not be comparable across papers.
- Automatically extracted numeric values require manual review.
- Culture-collection availability changes over time and must be manually confirmed.

## Modeling limitations

- CatPred-style predictions are enzyme-level priors, not strain-level production guarantees.
- Dextran product quality depends on enzyme kinetics, enzyme expression/secretion, medium, pH, temperature, substrate, process trajectory, and downstream processing.
- Kinetic or Simscape-style modeling requires time-series data; endpoint-only data support only black-box or reduced surrogate models.

## Open-science limitations

- Raw data from suppliers, companies, culture collections, or private NMR/GPC reports may not be shareable.
- Full protein sequences and raw API dumps should be reviewed before being included in public artifacts.
- The repository should separate reproducible code and schemas from private or licensing-sensitive data.

## Decision limitations

AESSP outputs are intended to choose pilot-screening candidates. They should not be used as final strain, enzyme, or process recommendations without experimental validation.
