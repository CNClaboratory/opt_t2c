# opt_t2c : Optimal type 2 criterion setting in Signal Detection Theory

## Summary

This is a toolbox for analyzing optimal type 2 criterion setting in classical Signal Detection Theory (SDT), as well as extensions to SDT that contain additional computational process influencing metacognitive decision making. It considers optimal type 2 criterion setting under several different optimization objectives: maximizing type 2 accuracy, maximizing type 2 reward, calibrating confidence to a threshold level of accuracy, and maximizing the difference between type 2 hit rate and type 2 false alarm rate.

This toolbox is a companion piece to a manuscript by Brian Maniscalco, Lucie Charles, & Megan Peters titled "Optimal Metacognitive Decision Strategies in Signal Detection Theory" (in review). The code can be used to reproduce and alter the analyses in that manuscript, as well as to conduct novel analyses.

## Description of files and folders

### Main functions

The root folder of the toolbox contains the main functions used to conduct the analyses.

- `opt_t2c` is a function that uses equations derived from classical SDT to compute optimal type 2 criteria, given an optimization objective and parameter settings.
- `opt_t2c_sim` is a function that uses computational simulations to find optimal type 2 criteria in models that extend classical SDT by allowing for type 2 evidence samples to differ from their type 1 counterparts. The simulation can be performed for any such possible model, provided that a function describing the model is provided. Details of how to do this are provided in the `opt_t2c_sim`'s documentation, and examples of implementing this functionality are provided in the simulation scripts described below.
- `t2noisySignalLoss` is a function that implements one such extension to classical SDT. In this model, type 2 evidence samples can undergo signal loss and receive additional noise. (Note that this function can implement both the "type 2 noise" and "type 2 signal loss" models described in the manuscript, hence the combined name "noisy signal loss".)

### Analysis scripts

The `analysis` folder contains scripts that use `opt_t2c` to plot the effects of various parameters (such as d' and c) on optimal type 2 criterion setting under various objectives.

- `opt_t2c_t2acc` investigates the effect of d', c, and p(S2) on the optimal type 2 criteria under the objective of maximizing type 2 accuracy. This analysis mirrors Figure 4 in the main manuscript.
- `opt_t2c_t2reward` investigates the effect of d', c, and Q2 on the optimal type 2 criteria under the objective of maximizing type 2 reward. This analysis mirrors Figure 5 in the main manuscript.
- `opt_t2c_calibration` investigates the effect of d', c, and p(correct)T on the optimal type 2 criteria under the objective of calibrating confidence to a threshold level of accuracy. This analysis mirrors Figure 6 in the main manuscript.
- `opt_t2c_hr2_far2` investigates the effect of d' and c on the optimal type 2 criteria under the objective of maximizing HR2 - FAR2. This analysis mirrors Figure 7 in the main manuscript.

### Simulation scripts

The `simulation` folder contains scripts that use `opt_t2c_sim` to plot the effects of various parameters (such as d' and c) on optimal type 2 criterion setting under various objectives when type 2 evidence samples differ from their type 1 counterparts. These simulations all use the "type 2 noisy signal loss" model, but can be readily adapted to investigate the behavior of other models as detailed in the documentation of `opt_t2c_sim`.

- `opt_t2c_sim_run` provides an example of running a single simulation using `opt_t2c_sim`.
- `opt_t2c_sim_sd2s` investigates the effect of type 2 noise on optimal type 2 criterion setting under various objectives. This simulation mirrors Figure 9 in the main manuscript.
- `opt_t2c_sim_ks` investigates the effect of type 2 signal loss on optimal type 2 criterion setting under various objectives. This simulation mirrors Figure 10 in the main manuscript.
- `opt_t2c_sim_t2reward` investigates how the effect of d', c, and Q2 on the optimal type 2 criteria for maximizing type 2 reward are modulated by type 2 noise and signal loss. This analysis complements and expands upon the analysis of Figure 5 in the main manuscript.
- `opt_t2c_sim_calibration` investigates how the effect of d', c, and p(correct)T on the optimal type 2 criteria for calibration are modulated by type 2 noise and signal loss. This analysis complements and expands upon the analysis of Figure 6 in the main manuscript.

## Notes

### Calibration under type 2 noise

To calibrate confidence to accuracy, the observer must find the location on the decision axis where estimated p(correct) equals the threshold level of accuracy, p(correct)T. It is ideal for the observer to use type 1 evidence samples for this purpose, but this has the consequence that type 2 noise does not affect calibration performance. Thus in `opt_t2c_sim` it is assumed that the observer uses type 2 evidence samples to estimate p(correct), in order to conduct an analysis on calibration under type 2 noise that is not trivial. Thus the results of these simulations reflect optimal type 2 criterion setting *given that* the observer estimates p(correct) using type 2 evidence samples.

In the plots produced by `opt_t2c_sim_sd2s` and `opt_t2c_sim_ks`, the plots for "estimated p(correct) at optimal c2" when "ignoring type 2 model (calc)" correspond to the observer's estimatation of p(correct) (as computed from type 2 evidence samples) at the location of the decision axis corresponding to the optimal type 2 criterion under classical SDT (which has no type 2 noise).

### Preserved qualitative patterns of optimal type 2 criterion setting under type 2 noise

In section 3.2.4 of the manuscript, we write:

> Complications arising from meta-d’ ≠ d’ might also be relaxed in cases where the researcher is more interested in investigating qualitative patterns in type 2 criterion setting, rather than in computing exact type 2 criterion values. For instance, it is likely the case that for most plausible models allowing for meta-d’ ≠ d’, the qualitative patterns that type 2 criteria should become more liberal with increasing d’, or more conservative with increasing Q2 or OT, hold. (Indeed, simulations confirm that these qualitative  patterns hold under type 2 noise. We do not show these data here, but provide code for running this analysis in the online code.) If the researcher is mainly interested in assessing such qualitative patterns, then it may not be necessary to engage with the complexities of characterizing exact type 2 criterion values when meta-d’ ≠ d’.

The analyses alluded to in this text are conducted in the scripts `opt_t2c_sim_t2reward` and `opt_t2c_sim_calibration`.
