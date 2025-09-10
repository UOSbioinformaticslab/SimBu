# Security Policy
Thank you for helping keep SimBu and its users safe.
This project is research/analysis software and does not run a hosted service, but we still take supply-chain and data-handling security seriously.

## 🔁 Supported Versions
Only the latest version is fully supported

| Version             | Supported          | notes              |
| -----------------   | ------------------ | ------------------ |
| vNext (main)        | :white_check_mark: | Security fixes land here first |
| Latest release      | :white_check_mark: | Receives security patches |
| Prior               | :x:                | Please upgrade     |


## Reporting a Vulnerability

Please include, when possible:
-	A clear description of the issue and why it’s a security risk.
-	A minimal reproducible example or PoC.
-	Affected versions/commits and your environment (OS, R version, package versions).
-	Any suggested mitigations or patches.

We aim to release fixes within 5 days when feasible.

## 🧭 Scope

In scope:
	•	Source code in this repository (R code, scripts, package metadata).
	•	Provided example data and configuration files.
	•	Build/CI configuration that impacts consumers of the package.

Out of scope:
	•	Vulnerabilities in third-party dependencies not patched upstream.
	•	Issues only reproducible with non-default, unsafe flags or modified code.
	•	Non-security bugs (please use Issues for those).

If you’re unsure whether something is in scope, report it privately and we’ll triage.

## 🔁 Supported versions
We generally support the most recent release and the previous minor version.

## 🧩 Dependencies & supply chain
	•	Dependency updates are tracked via GitHub’s Dependabot and advisories.
	•	We prefer upstream fixes for third-party packages; where necessary we may apply temporary pins/workarounds.
	•	Release artifacts are built from tagged commits; avoid installing from untrusted forks.

If you identify a vulnerable dependency path (e.g., via renv::snapshot() or pak), please include the lockfile/sessionInfo() and the CVE/CVSS reference.

## 🧪 Reproducible security reports

When reporting, please attach a minimal script or R Markdown showing:
```r
# Example scaffold
sessionInfo()
# install.packages("SimBu") or devtools::install_github("stef1949/SimBu")
library(SimBu)
# minimal code that triggers the issue
```
