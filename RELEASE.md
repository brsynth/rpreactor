# Release history

## 0.7.1
- docs: use nbsphinx to showcase usage
- feat: silence RDKit warnings
- fix: sqlite3 IntegrityError with new results and commit
- fix: remove *all* rule_mol already yield in compute()

## 0.7.0
- docs: sphinx documentation is available
- fix: processes don't hang anymore!
- refactor!: replace multiprocessing by pebble
- refactor!: JSON export in the CLI + compute as a generator
- refactor!: compute couples of rule x chemicals
- feat!: sqlite3 as a cache database + return rdkit Mol object
- feat!: create sqlite3 database from retrorules
- feat: small optimization of memory usage
- feat: support for stoichiometry
- feat: compute on subsets of rules or chemicals  
- feat: fusion project rpchemtools
- ci: package publication is automated by github workflows

## 0.6.2
- docs: install instructions update
- chore: introduce RELEASE file

## 0.6.0
- BREAKING CHANGE: move core methods to Core script
- docs: fix set up development environment

## 0.5.0
- fix: clean-up Pool after compute()
- chore: CI with github action
- doc: updated README


