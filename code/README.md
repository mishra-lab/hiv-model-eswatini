TODO: needs major update

# Notation Conventions
- dimensions: see __init__.py
- variables:
  - X: modelled population
  - P: parameter set (dict)
  - T: targets (list)
  - R: model result (dict)
  - Ps,Ts,Rs: list of any above
- parameter names:
  - *_{stuff}: stuff is either the dimensions that * is stratified by (e.g. si: sex and activity)
               or a specific value within * (e.g. FSW_2010: female sex workers in 2010)
  - R*: relative value of *
  - dx,tx,vx: rates of diagnosis, treatment, and viral suppression
- _ = None, a shorthand for adding new dimensions to ndarrays

# Structure
- model - mostly everything related to running the model
  - __init__: common stuff for slicing by name & "global" definitions of dimension stuff
  - check:    [old] model checking stuff
  - foi:      subset of system stuff just related to force of infection (incidence)
  - fit:      easily create many plots for debugging model parameters
  - main:     where everything starts
  - out:      main & supporting functions for calculating outputs from the model
  - params:   functions for defining & sampling parameters
  - plot:     functions for plotting stuff, often calling out functions to compute it too
  - scenario: working functions to run the main-ish analyses
  - system:   all non-foi functions related to the ODE system, plus some helpers for solving
  - target:   definitions of calibration target & infrastructure to compare with model outputs
              n.b.: some outputs/targets are based on one "pop"ulation, some based on two
              the "vs" stuff helps organize two-population outputs, i.e how to compare them ("vsop")
- utils - python and R "utility" (helper) functions
  - ops: a lot of tiny helper functions
  - stats: some wrappers for scipy.stats 
  - tarray: custom class for time-varying parameters, including ndarrays
    - arguments: time points (list), known values or "knots" (ndarray), with t dimension last
    - implementation: we fit a monotonic spline for every array value (separately) at initialization
    - the fitted spline parameters are stored, and then can be evaluated as needed
    - use A(t) to evaluate all the splines at "t", then reshape the list to the appropriate shape
  - deco: we use decorators to evaluate standard code before/after a bunch of functions
    - e.g. rmap selects which keys within "R" should passed as keyword arguments to the function
    - e.g. tslice slices some time-varying arguments at "t", assuming those args have t-dim first
    - e.g. nanzero sets all nans in the return argument to zero
- params - upstream code for parameterization
  - distr.py: helper code for estimating distribution parameters from mean +/- 95% CI
  - fsw: main analysis code for FSW risk stuff

# Platform
everything built & tested on:

- linux mint 20.1
- python 3.8.10
  - we depend on ordered dicts (3.7+)
- R 3.6.3
- atom 1.57.0

... no guarantees things will run elsewhere

