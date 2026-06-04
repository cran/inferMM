# inferMM 0.0.3

- Removed the `nls` fitting path so all single-curve fits now use the manuscript's Michaelis-Menten profile-score estimator.
- Updated clustered initialization to use the pooled single-curve profile-score fit without an `nls` fallback.
- Refreshed package documentation to reflect the unified profile-score implementation.

# inferMM 0.0.1

- Initial development version.
- Added variance-aware Michaelis-Menten fitting and screening functions.
- Added grouped workflows, simulation helpers, tests, and demo data.
