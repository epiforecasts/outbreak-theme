# Project-Specific Claude Instructions

This project follows the Abbott et al. "Workflow for Infectious Disease Modelling" (11 steps, 00-10).
Reference: https://github.com/seabbs/infectious-disease-modelling-with-multiple-datasources

## Workflow Step PRs

- Each workflow step (step00, step01, ..., step10) gets its own branch targeting main
- Each PR should ONLY edit its own step file (e.g., step02 branch only edits step02_*.md)
- Branch naming: `step00`, `step01`, etc.

## CodeRabbit Reviews

- CodeRabbit automatically reviews PRs
- Address actionable comments; nitpicks are optional
- Reply to inline comments in the PR thread after fixing (use `gh api .../comments/{id}/replies`)
- CodeRabbit should auto-approve once concerns are addressed

## PR Descriptions

- NO test plan section
- Keep descriptions concise: summary bullets only
- Link to issues if applicable

## Julia Code

- Model implementation in `hpai-challenge/` directory
- Use Turing.jl for probabilistic programming
