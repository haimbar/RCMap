# RCMap 0.3.1

## Improvements

- **Cluster names in statement report.** `StatementSummaryNN.csv` now includes
  a `ClusterName` column (last column) containing the user-assigned cluster
  name for each statement row, in addition to the existing `ClusterNo` integer
  column.

- **Cluster names in ANOVA and Tukey output.** The one-way ANOVA table
  (`ANOVANN.txt`) and Tukey HSD post-hoc tests (`TukeyNN.txt`) now label the
  cluster factor levels with user-assigned cluster names rather than integer
  indices. Pairwise Tukey comparisons therefore appear as e.g.
  `"Health-Education"` instead of `"2-1"`. When no names have been assigned the
  cluster index is used as before.

# RCMap 0.3.0

## New features

- **Fuzzy pile label matching.** Pile labels entered by different sorters are
  now automatically grouped using Jaro-Winkler string distance, so typographical
  variants such as "Health", "helath", and "hlth" are merged into a single
  canonical label before cluster name suggestions are computed. The matching
  threshold is configurable via `config.txt` (`fuzzy_label_threshold`, default
  0.15).

- **Interactive canonical label editor.** A new Settings menu option (option 5,
  "Edit pile label canonical names") lets users review all canonical labels
  alongside the original sorter-supplied labels. Labels that were automatically
  merged are highlighted in yellow. Any canonical label can be renamed
  interactively; changes are persisted to `dictionaries.RData` and
  `CMapSession.RData` immediately.

- **Project configuration file (`config.txt`).** On first load RCMap creates a
  plain-text `config.txt` in the project folder listing every tunable parameter
  with its default value and a short description. Missing settings are appended
  on subsequent loads; values already present are never overwritten. Supported
  settings: `ratingscale`, `clust_method`, `dist_metric`, `color_scheme`,
  `n_clusters`, `mds_seed`, `splithalf_seed`, `splithalf_B`,
  `fuzzy_label_threshold`, `jaccard_threshold`.

- **Configurable seeds and thresholds.** The MDS jitter seed (`mds_seed`), the
  split-half seed and replication count (`splithalf_seed`, `splithalf_B`), and
  the Jaccard instability threshold (`jaccard_threshold`) are now read from
  `config.txt` and stored in the session, making all analyses fully
  reproducible without code changes.

- **Non-consecutive IDs supported.** Statement IDs, sorter IDs, and rater IDs
  no longer need to be consecutive integers; any unique numeric identifiers are
  accepted.

- **Graceful rater ID mismatch handling.** If a rater appears in the ratings
  file but not in the demographics file, a yellow warning is shown and analysis
  continues. If a rater appears in the demographics file but has no rating data,
  they are automatically excluded with a warning.

- **Error recovery in folder selection.** If the selected data folder is missing
  required files, an error message is displayed and the program returns to the
  main menu instead of terminating.

## Bug fixes

- Settings menu no longer became unresponsive after returning from a submenu
  (stale `showSettingMenu` guard caused the submenu to be silently skipped).
- Weight loop now correctly matches sorters by their original ID rather than
  their sequential index, fixing incorrect weighting when sorter IDs are
  non-consecutive.
- `nCards` in `getAdjMatrices` is now derived from the statement count rather
  than the number of distinct card values seen in pile data, preventing
  out-of-bounds writes when some statements were never sorted.
- Rating range validation in `read_ratings` now correctly coerces values to
  numeric before comparison (previously used character matrix comparison, which
  silently passed invalid values).
- `dataDir` NA check in `initCMap` now guards against `nzchar(NA)` returning
  `NA` rather than `FALSE` when the user cancels the folder picker.
- `Weights.csv` is now fully validated (column presence, row count, numeric
  values, positive values) before use.
- Extra raters in `Demographics.csv` (those with no corresponding rating rows)
  are now actually filtered out rather than only warned about.
- `sorterReport` now prints original sorter IDs instead of sequential indices.
- `class() == "factor"` comparisons replaced with `inherits()` throughout.

## Documentation

- Full roxygen2 documentation added for all previously undocumented functions:
  `read_piles`, `read_ratings`, `read_demographics`, `initCMap`, `get_cmat`,
  `clusterings`, `toText`, `colLab`, `topLine`, `rcmenu`, `clusterConfig`,
  `clusterNames`, `update_clusters`, `loadRCMapData`, `choose_dir`,
  `choose_directory`, `getOS`, `RCMapMenu`, `plotjaccard`.
- Several function descriptions corrected (e.g. `splitHalf` had copy-pasted
  description from `MPindex`; `plotjaccard` described a "MADD dotchart").
- `Documentation.Rmd` updated with a new "What's New" section, expanded input
  data guidance, updated Settings menu, and a new "Project Configuration File"
  section documenting all `config.txt` settings.
- Added `stringdist` to package dependencies.

# RCMap 0.2.3

- Bug fixes and stability improvements.
- Sorters and raters no longer required to be the same group.

# RCMap 0.2.0

- Major rewrite: input data now expected as four plain-text CSV files
  (`Statements.csv`, `SortedCards.csv`, `Ratings.csv`, `Demographics.csv`).
- New menu-driven command-line interface (`RCMapMenu()`).
- Split-half reliability analysis added.
- Misplacement index added.
- Parallel coordinates plot (formerly "pattern matching") supports more than two
  cohort/variable combinations.
- Binary cohort variables automatically derived from quantitative demographics
  (X > median(X)).
- Choice of Euclidean or hyperbolic distance metric for MDS and clustering.
- Two automatic cluster-count methods: within-cluster sums of squares and
  average silhouette (via `factoextra::fviz_nbclust`).
- Dot-chart and bar-chart rating plots added.
- Colored dendrograms and phylogenic tree visualisation added.
- Two color schemes: `rcmap` (up to 21 clusters) and `rainbow`.
