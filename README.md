# skill_game

This repository holds the analysis functions for the video gaming skill study. 

Analyses Content:
LEARNING:
(1) Distance Travelled
(2) Hazard
(3) Kinematic variability
(4) Policy deviation

PROBE:
(1) Distance Travelled
(2) Hazard
(3) Kinematic variability
(4) Policy deviation
(5) Interaction (kinematics vs. success)
	(a) Policy deviation & Distance Travelled

----- PROBE(5)(a) There is a concern in splitting the analysis of policy 
deviation into successful runs and failure runs. When considered 
together, the overall average policy deviation worsens substantially 
during probes. This was almost innevitable given that the distance 
travelled dropped by 10%. But analyzing deviation separately for 
successes and failures seems to only introduce a selection bias. If we 
only look at trials that were good, then it's of no surprise that we 
would see good kinematics. (N.B.: this is an issue with kinematic 
variability as well, since we are only measuring kinematic variability 
from trials that succeeded. This, of course, is necessary because we 
cannot measure the variability of data we don't have). At least two 
things make this analysis interesting, and indeed valid. The first is 
that there is dynamic range in both the failure trials and the 
successful trials across learning. Thus, it is not as if successful 
trials invariably have the same lower level of policy deviation - they, 
too, show signs of learning, and same with failure trials. In other 
words, there is room for these trials to deteriorate to, namely, a 
previous state of learning. This is especially true of the failure 
trials, which show greater change across learning. And despite this, 
policy deviation does not change during probes, given the status as 
failure or success.

Nevertheless, to better demonstrate this point, two additional analysis 
will be done. The first is to show that the way in which failure trials 
fail is that they start out looking like successes, but then take a 
sharp deviation leading to their demise. This can be done by showing the 
policy deviation along the length of the track for both successes and 
failures, and showing how the signal for failures increases as they 
trail off.

A second way to do this is to compute the average difference in distance 
travelled between trials in the probe and pre-probe windows who have 
matched kinematics up to a certain point. The point used to determine 
that kinematics are matched can be changed across the track so that the 
difference in distance travelled can be computed all along the track. 
Probe trials should have lower distance travelled than pre-probe trials 
despite having the same kinematics up to a given point.This analysis 
shares the same spirit as a hazard analysis. It also would further 
confirm the point from the above analysis that failure trials are 
actually successes until they are not.

