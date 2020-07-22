# rpreactor

Apply reaction rules and parse results

## Installation
```bash
# Setting up the conda environment
conda create --name myenv python=3
source activate myenv
conda install --channel rdkit rdkit=2019.03.1.0
conda install --channel tduigou rpchemtools
```

## Usage

From command line:
```
python -m rpreactor.cli --with_hs true inline --inchi "InChI=1/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)" --rsmarts "([#8&v2:1](-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11])-[#1&v1:12])>>([#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])(-[#8&v2]-[#1&v1])-[#8&v2:1]-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11].[#7&v3](=[#6&v4]1:[#7&v3]:[#6&v4](-[#8&v2]-[#1&v1]):[#6&v4]2:[#7&v3]:[#6&v4](-[#1&v1]):[#7&v3](-[#6&v4]3(-[#1&v1])-[#8&v2]-[#6&v4](-[#6&v4](-[#8&v2]-[#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1:12])(=[#8&v2])-[#8&v2]-[#1&v1])(-[#1&v1])-[#1&v1])(-[#1&v1])-[#6&v4](-[#8&v2]-[#1&v1])(-[#1&v1])-[#6&v4]-3(-[#8&v2]-[#1&v1])-[#1&v1]):[#6&v4]:2:[#7&v3]:1-[#1&v1])-[#1&v1])"
```

## Build
```bash
conda create --name build_env python=3
conda install conda-build
conda build recipe --channel rdkit --channel tduigou
```

## TO DO

- Add code to remove stereo in Utils.standardize_chemical, update tests accordingly.

## Buggy tales

### 20190110.01 -- Timeout issue
  * Checking results: timeout neither reached... Maybe fine but strange.
  * Checking results: bunch of error messages mentioning unknown "pool" variable (e.g. `"match_error": "name 'pool' is not defined"`)
  * Bug: typo led to wrongly handling Pool of process
  * Fix: commit ca9ff17f1e @ rule_fire
  * Imply: rerun jobs impacted (DONE, all jobs have been reran)

### 20190110.02 -- Empty InChI for products
  * Checking results: some inchi string are empty, e.g.: rule (radius 0) `8a4666fc9dab` applied on `MNXM366`
    - `8a4666fc9dab`: `([#6&v4:1]1:[#6&v4:2]:[#6&v4:3]:[#6&v4:4]:[#7&v3:5]:[#6&v4:6]:1)>>([#6&v4:1](-[#6&v4:2]-[#6&v4:3]:[#6&v4:4]-[#7&v3:5]-[#1&v1])(-[#6&v4:6](-[#7&v3](-[#1&v1])-[#1&v1])-[#1&v1])-[#1&v1])`
    - `MNXM366`: `[H][O][c]1[c]([C]([H])([H])[H])[n][c]([H])[c]([C]([H])([H])[O][P](=[O])([O][H])[O][H])[c]1[C]([H])([H])[N]([H])[H]`
  * Reason: sanitization issues of products (see below) leading to empty InChI.
  ```
  [16:06:30] non-ring atom 2 marked aromatic
  10/01/2019 16:06:30 -- WARNING -- Partial sanization only
  [16:06:30] ERROR: Unrecognized bond type: 0
  [16:06:30] non-ring atom 2 marked aromatic
  10/01/2019 16:06:30 -- WARNING -- Partial sanization only
  [16:06:30] ERROR: Unrecognized bond type: 0
  ```
  * Reason: probably due to rule miss-encoding. Below are rule known to raise with bug (checked for radius 1, 3, 5, 7, 9, 10):
  ```
  Count Rule_hash
  15 24f4a1100077
  57 300bddc687b7
  15 88eee4c1d443
  57 e9f88e91236e
  ```
  * Consequences: probable wrong transformations.
  * Workaround: reject all results from a couple substrate-rule as soon as an "empty" InChI is produced
  * Fix: commit e980692922 @ rule_fire
  * Imply: take care of this situation in results already generated / rerun code and make new results (DONE: new results generated)

### 20190110.03 -- Cannot convert back InChI to RDKit mol
  * Rule (radius 0) `3012ddcc4356` applied on `MNXM2183`
    - `3012ddcc4356`: `([#6&v4:1](-[#7&v3:2]:[#6&v4:3](:[#7&v3:4]):[#6&v4:5])-[#8&v2:6])>>([#6&v4:1]=[#8&v2].[#7&v3:2]:[#6&v4:3](:[#7&v3:4]-[#1&v1])-[#6&v4:5].[#8&v2:6]-[#1&v1])`
    - `NXM2183`: `[H][N]=[c]1[n][c]([O][H])[c]2[n][c]([H])[n]([C]3([H])[O][C]([H])([C]([H])([H])[O][H])[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[c]2[n]1[H]`
  * Reason: issue with valencies, aromatic cycle and nitrogen.
  * Fix: none at the moment, but should probably be handled at the result parsing step (i.e. after rule firing)

### 20190112.01 -- The unkillable Pool
  * Execution hang (sometimes) when renewing the Pool (done after every time out)
  * Due to timed out job lock that cannot be acquired (deadlock)
  * Issue do not happen every time
  * 100% chance of hang when the time out == 0. Hypothesis: the job was not fully initialized, cannot acquire a lock that have not been set up
  * But hang can also happen if more "reasonable" timeout condition (e.g. 5 seconds). Hypothesis: the job have been killed while it was actually ending, something went wrong.
  * Straightforward fix for all cases: force to release the job lock each time. Commit 52466a536e @ rule_file
