from model.strat import strats

tol = 1e-8

# dimension reference
# (X) Population:
# - (s) Sex:               0: Women, 1: Men
# - (i) Activity:          0: Lowest Activity, 1: Medium Activity, 2: Lower Risk Sex Work, 3: Higher Risk Sex Work
# - (k) Inert Partnership: 0: None, 1: Main/Spousal, 2: Casual, 3: One-Off Sex Work, 4: Repeat Sex Work
# - (h) Health:            0: Susceptible, 1: Acute HIV, 2: CD4 > 500, 3: 350 < CD4 < 500, 4: 200 < CD4 < 350, 5: CD4 < 200 (AIDS)
# - (c) Care:              0: Undiagnosed, 1: Diagnosed, 2: Virally Unsuppressed, 3: On ART, 4: Virally Suppressed
# (*) Other:
# - (a): Sex Act Type:     0: Vaginal, 1: Anal
# - (p): Partnership:      0: Main/Spousal, 1: Casual, 2: One-Off Sex Work, 3: Repeat Sex Work
