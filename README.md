# dr1-wl-testsuite
A repo for the collection of routines to carry out tests on the shear catalogues for DR1. Contains both WL1 and WL2 contributions.

# Usage

Place any routines necessary to produce plots in the `modules` directories: this can include anything that returns a figure / axes, as well as modules which produce necessary data for plots. Place scripts to output plots in the `plotters` directory. If possible, please try to use a naming convention that ties to specific tests, as detailed in the WL-1 and WL-2 test document.

# Test definitions

The testing strategy is given [here](https://docs.google.com/document/d/1CjEPUxUHbXvhBhjemsJGFzE79WinPhtVvm62RUHN16E/edit?tab=t.0#heading=h.zahu9s6fxywk).

See [TR1/DR1 test lists (redmine)](https://euclid.roe.ac.uk/projects/dr1-kp-jc-2/wiki/List_of_tests_for_DR1_data_products). The [WL test document](https://docs.google.com/document/d/1-NlqlUvYEPviI4lITrRwDtV6cryN1ydo63bwU_zhJto/edit?tab=t.0#heading=h.2w9qeojvo2op) is supplementary: this can also include tests that don't **need** to be run on each dataset / cycle.
