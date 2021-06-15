%% Chronic carriers. How many did Ames' 1943 study miss.
% In that study, an individual was a chronic carrier if they could not
% produce two consecutive specimens 

% What's the probability that a chronic carrier doesn't fulfill the release
% criteria in a year? Each time, they can produce an N or P if the last
% test was P, but they must produce a P if the last one was N. For the
% moment, let's assume that they got a test every single week. (It may
% matter how many tests they actually tried, no?)

% According to Gilman 1975 study, a rectal swab has a sensitivity of 23/62
% or 21/60 if you want to think of the bone-marrow-positive patients only.

% So the sensitivity parameter is betarnd(1, 23+1, 62-23+1)
% Get N_iterations of a 49-element vector of Ns and Ps. Bernoulli random
% process with a parameter given by the beta distribution above.
% For each 49-element vector, let's see if the chronic carrier "got away"


% The reason I think this is important is that it could affect our estimate
% of the contribution of chronic carriers, and it could affect the dynamics
% of the model a bit because a larger proportion of the adults would be
% "stolen" away in the C compartment where they cannot get sick. The
% dynamics could work differently.

% Let's leave this re-analysis for the dose-response/identifiability paper.
% Then we can talk about the implications of all the unknowns about chronic
% carriers (as a secondary goal of the paper).