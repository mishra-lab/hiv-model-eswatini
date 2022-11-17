(5 minutes)

Thanks for checking out our posted titled,
*Beyond instantaneous partnerships: capturing partnership-level herd effects in compartmental models of sexually transmitted infections*

**Background**
- compartmental models of sexually transmitted infections, including HIV, are key tools in epidemic response
- such models define the rate of new infections using the average infection prevalence in the modelled population
- so, newly infected individuals immediately contribute to more infections
- in other words, we have an assumption that partnerships are effectively instantaneous
- in reality, partnerships include periods before and after transmission
- previous work has shown that this hidden assumption can bias estimates of intervention impact

**Objectives**
- we sought to develop a new compartmental model which avoids this instantaneous partnership assumption
- then, we wanted to compare the modelled incidence under the old instantaneous model, versus the our new proposed model:
  specifically incidence overall, and incidence within longer vs shorter partnerships

**Instantaneous Partnerships**
- before developing our new approach, let's look closer at how instantaneous partnerships are actually modelled
- the incidence rate includes a number of terms, defined here, including the prevalence of infection,
  and in fact the partnership duration
- so how is this duration used?
- in the real world, transmission may occur within a partnership after a number of sexual contacts
- however, the partnership might then continue for days, months, or even years after transmission
- all subsequent contacts within the partnership cannot contribute to transmission,
  so these contacts are sometimes called "wasted" contacts
- so, within this incidence equation, the average rate of transmission for each partnership
  is reduced to account for these wasted contacts
- however, we're forced to make one of two assumptions:
  1) if we use the true partnership duration, we underestimate early transmission by spreading out that risk over a longer period,
- on the other hand,
  2) if we assume a maximum duration of 1-year, which is commonly done, then we underestimate the impact of wasted contacts.
- so the two limitations of this approach are:
  1) that partnerships are instantaneous, and
  2) that we're forced to capture either early transmission or wasted contacts, but we cannot capture both

**Proposed Model**
- so how can we develop a model which properly accounts for partnership duration? (Objective 1)
- well, the core idea is to recognize that individuals who recently acquired or transmitted infection
  are likely still with the same partner
- so, they cannot contribute to transmission
- to keep track of these individuals, we introduce a new compartment, or "state" (shown in yellow)
- then, we remove these individuals from the incidence rate equation,
  until they change partners, as defined by the inverse partnership duration
- it gets a little more complicated when individuals have multiple concurrent partnerships, or different partnership types
- in this case, these individuals are still included in the incidence rate equation,
  but their effective numbers of partners are reduced by one
- and this proposed model represents the solution to our objective 1

**Simulation Study**
- Next, we integrated the proposed model into an existing model of heterosexual HIV transmission in eSwatini
- the model includes over 200 compartments, including 8 different risk groups, 5 stages of HIV, and 5 states related to treatment
- the model also includes 4 partnership types, including longer "main" partnerships and shorter "casual" and "sex work" partnerships
- We calibrated this model to eSwatini data on HIV prevalence, incidence, and treatment,
  first under the instantaneous model, and again under the proposed model (results in figure)
- then, for objective 2, we conducted two experiments:
  (a) we compared overall incidence under the two models using the same parameters
  - this allows us to see the direct impact of the instantaneous partnership assumption
  (b) we examined incidence within longer and shorter partnerships, this time using model-specific parameters (from the two calibrations)
  - this allows us to see the indirect impact on which types of partnerships are prioritized for prevention

**Simulation Results**
- here we can see the results of experiment (a)
- with equal parameters, overall incidence under the old instantaneous model is substantially higher than under the proposed model
- these differences grow over time,
  reflecting the accumulation of infections "trapped" in longer partnerships after transmission in the proposed model,
  or what we call "partnership-level herd effects"
- then, in experiment (b), we see that much more transmission occurs via longer partnerships
  under the instantaneous model versus under the proposed model
- and conversely, less transmission occurs via shorter partnerships

**Implications**
- in summary, our proposed model overcomes a decades-old partnership modelling problem
  - we can use true partnership duration and thereby capture both early transmission and wasted contact dynamics
  - we also capture the accumulation of partnership-level herd effets
  - in doing so, we avoid the need for more complex modelling frameworks, such as network models.
- we also have shown how existing models of HIV and other STI likely
  overestimate the impact of prevention in long-term partnerships and
  underestimate the impact of prevention in short-term partnerships

**References**
- all of our code is freely available on github. Thank you for your attention


