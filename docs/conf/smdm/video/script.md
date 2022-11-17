(2 minutes)

Hello and thanks for checking out our poster, titled
*Beyond instantaneous partnerships: capturing partnership-level herd effects*
*in compartmental models of sexually transmitted infections*

- Compartmental models of infectious disease are *key tools* in epidemic response,
including for sexually transmitted infections like HIV.
- These models define the rate of new infections using the *average infection prevalence*
- So, *newly infected individuals immediately contribute to more infections*
- This effectively assumes that sexual partnerships are instantaneous,
and previous work has show that this hidden assumption can *bias estimates of intervention impact*

So, we sought to develop a new compartmental model which avoids this instantaneous partnerships assumption.

- The core idea of our model is to recognize that individuals who *recently acquired or transmitted infection*
might be still with the same partner, and so they *might not contribute to incidence*.
- To keep track of these individuals, we introduce a new compartment,
and we remove them from the incidence rate equation until they form a new partnership.
- This approach can be further generalized for multiple partnerships.

- We integrated our new model into an existing model of heterosexual HIV transmission in eSwatini,
and compared simulated epidemics under the old instantaneous model versus our new model.
- We found that the old instantaneous model overestimates incidence, because it fails to capture
when infections become "trapped" in longer partnerships after transmission,
or what we call *partnership-level herd effects*.
- The old instantaneous model also *overestimates the impacts of prevention in longer partnerships*
and underestimates the impacts of prevention in shorter partnerships.

- Such findings have notable implications for
the existing body of model-based evidence for STI and HIV epidemic response.
- We hope that our new model can be incorporated into transmission models going forward,
and you can find our code at the github link show here.

Thank you,