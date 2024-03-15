# Eswatini HIV Transmission Model: Code

## overview

- model:    all model implementation & scenarios code (py)
- params:   ad-hoc analyses to inform model parameters (R & py)
- post:     ad-hoc analyses on model outputs, mainly stats & plots (R)
- profile:  enable profiling & run model 10 times without parallel (py)
- tikz:     diagram-drawing code (latex/tikz)
- toy:      ad-hoc analyses of toy systems (py)
- utils:    collection of utility functions (R & py)
- makefile: commands for running all major analyses

## model code

- __init__: global model config stuff
- debug:    minimal code to run the model & plot some figures
- fit:      functions for plotting model outputs vs calibration targets (does not run fitting)
- foi:      force of infection functions, including transmission-prob, mixing, & incidence
- out:      functions for computing (stratified) model outputs
- params:   functions to specify & sample model parameters (inputs)
- plot:     functions for plotting model outputs & calibration targets
- strat:    a collection of strata meta-data for plotting in python
- system:   main functions for running the model
- target:   functions to specify & check model outputs vs calibration targets
- tpaf:     function to compute TPAFs by running model with vs without masked transmission
- scenario: files for running the model for main analyses
  - __init__: global config stuff, including model version & filenames
  - main:     bash script for all major analysis steps & running on SciNet
  - imis:     implements adapted IMIS calibration
  - art:      runs recalibration, outputs, & sensitivity analyses for ART inequalities paper
  - foi:      runs recalibration & outputs for FOI equations paper

### model code notation

- X:    # people in each compartment over time (model state variable)
  - Xk: X without summing across EPA partnership dimension
- T:    list of calibration targets
  - Ti: individual calibration target
- P:    dict of model parameters (inputs)
  - Ps: list of P dicts, mainly for running model multiple times
- D:    dict of sampling distributions for P
- R:    dict of model result/return data
  - Rs: list of R dicts, mainly from running model multiple times
- key config objects:
  - N:     dict of IMIS sample counts (numbers of model runs)
  - b:     int index of this IMIS batch
  - tvec:  dict of time vectors for different model runs
  - uid:   ID for this version of model, params, & targets
  - nid:   ID for current N

### model code conventions

- dimensions & strata: see model.__init__ (comments)
- targets & outputs: see model.out (labels)
- parameter names (see also code/params/tab/par.defs.csv)
  - x_a:   parameter 'x' is either stratified by dimension 'a'
           or corresponds to a specific value 'a'
  - Px_a_b: proportion of 'x' which is 'b' among 'a'
  - Rx_a:b: relative value of 'x' among 'a' vs 'b'
  - aRx:   additional relative value of 'x', like Rx = 1+aRx
  - C12m:  reported partners in past 12 months
           not equal to change rate due to duration adjustment
  - K:     current number of partners (cross-sectional)
  - Q:     partnership change rate (per-year)
  - F:     frequency of sex per partnership (no dilution)
  - dur:   duration in risk group or partnership length
  - turn:  rate of turnover among risk groups
  - lpref: log-preference of mixing among risk groups
  - beta:  probability of transmission per sex act
  - gud:   genital ulcer disease (increases both sus & inf)
  - prog:  rate of HIV progression among stages
  - dx:    rate of HIV diagnosis among PLHIV
  - tx:    rate of ART initiation among diagnosed
  - vx:    rate of viral suppression among those on ART
  - unvx:  rate of ART failure/discontinuation among VLS
  - revx:  rate of re-viral suppression among fail/discontinued

## post code

- config:    global config stuff, including model version & filenames
- post:      functions to support loading & plotting posterior param distributions
- tpaf:      functions to support plotting TPAFs
- wiw:       functions to support plotting out.wiw (who infected whom)
- main.base: functions for analysis of base calibration
- main.art:  functions for plotting & analyses for ART inequalities paper
- main.foi:  functions for plotting & analyses for FOI equations paper

## utils code

- ops:    small helper functions (py & R)
- stats:  functions to support stats; mainly remap args for scipy.stats (py)
- fio:    functions to save & load data as: .npy, .csv, .txt (py)
- deco:   decorators to execute code before and/or after other functions (py)
- tarray: subclass of np.ndarray to interpolate time dimension on-the-fly (py)
- plot:   functions to support plotting, including from out.expo (R)

## tikz code

- makefile: lists all figure names (un-comment to run a figure)
- maketikz: shell script to compile main to .pdf, crop, & save in out/fig/tikz/
- main:     skeleton latex file to compile
- config:   packages & macros to support other figures
- colors:   standardized color definitions to match strata
- * :       (all other files) tikz code for each specific figure

## platform

everything built & tested on:

- linux mint 20.1
- python 3.8.10
- R 3.6.3

windows plebs may have trouble; maybe try WSL

## authorship

model developed & maintained by:

- Jesse Knight [2022.01 - 2024.03] jesse.x.knight@protonmail.com
- Siyi Wang [2024.04+]

see docs/ * /meta.tex for contributors & co-authors

