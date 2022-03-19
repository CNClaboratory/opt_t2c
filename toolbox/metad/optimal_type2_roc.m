function [f2 h2] = optimal_type2_roc(d, c, pS2, grain, bound)

% [f2 h2] = optimal_type2_roc(d, c, pS2, grain, bound)
% 
% Given the specifications of an SDT model, estimate the optimal overall 
% type 2 ROC curve based on type 2 likelihood ratio.
% 
% inputs
% ------
% d, c - d' and c on the equal variance SDT model
% pS2 - Prior probability of displaying an S2 stimulus.
%    If not specified, defaults to 0.5.
% 
% grain, bound - parameters controlling construction of the ROC curve.
% This function computes type 2 likelihood ratios at different samples of the 
%    type 1 decision axis. The parts of the decision axis that are sampled are 
%    c-bound : grain : c+bound. So setting grain to smaller values divides 
%    the decision axis more finely, and setting bound to larger values probes 
%    more parts of the decision axis.
%    If  not specified, grain defaults to .01 and bound defaults to 10.
% 
% outputs
% -------
% f2 and h2 are vectors containing overall type 2 false alarm rates and hit 
%    rates generated by using the type 2 likelihood ratio as the type 2 
%    decision axis. These are samples from the optimal overall type 2 ROC curve 
%    that can be constructed for the given d', c, and p(S2). They can be used 
%    e.g. to construct a plot of the type 2 ROC curve or to estimate area 
%    under the curve.


%% initialize parameters

if ~exist('pS2','var')
    pS2 = .5;
end

if ~exist('grain','var')
    grain = .01;
end

if ~exist('bound','var')
    bound = 10;
end


%% find response-conditional type 2 likelihoods

x_rS1   = c-bound : grain : c-grain;
t2L_rS1 = t2L(x_rS1, d, c, pS2);

x_rS2   = c+grain : grain : c+bound;
t2L_rS2 = t2L(x_rS2, d, c, pS2);

t2Ls = [t2L_rS1 t2L_rS2];


%% construct the overall type 2 ROC curve

f2 = zeros(size(t2Ls));
h2 = zeros(size(t2Ls));
for i = 1:length(t2Ls)
    
    %% find the response-conditional x values that generate matching type 2 likelihoods
    L = t2Ls(i);
    [xx ind_rS1] = min( abs(L-t2L_rS1) );
    [xx ind_rS2] = min( abs(L-t2L_rS2) );
    
    xL_rS1 = x_rS1(ind_rS1);
    xL_rS2 = x_rS2(ind_rS2);
    
    %% calculate overall type 2 F, H from matched-likelihood x values
    [f2(i) h2(i)] = sdt_2_ot2_roc(xL_rS1, xL_rS2, d, c, pS2);    

end

[f2 ind] = sort(f2);
h2       = h2(ind);


end


%% compute type 2 likelihood
% specifying d', c, and prob(S2 stimulus presentation)
% determines the type 2 likelihood ratio at x, LR2(x) = f(x|C) / f(x|I)
%
% equations follow Galvin et al (2003) table 4

function L = t2L(x, d, c, pS2)

pS1 = 1 - pS2;

S2mu =  d/2;
S1mu = -d/2;

fS1 = normpdf(x,S1mu,1);
fS2 = normpdf(x,S2mu,1);

pHit  = 1 - normcdf(c,S2mu,1);
pCR   = normcdf(c,S1mu,1);

pC = pS2*pHit + pS1*pCR;
pI = 1 - pC;

fC(x <= c) = fS1(x <= c) .* pS1 ./ pC;
fI(x <= c) = fS2(x <= c) .* pS2 ./ pI;
fC(x > c)  = fS2(x > c) .* pS2 ./ pC;
fI(x > c)  = fS1(x > c) .* pS1 ./ pI;

L = fC ./ fI;

end


%% compute overall type 2 (FAR, HR) from SDT model
% given a full SDT model, find the overall type 2 (FAR, HR)
% for a given pair of response-conditional type 2 criteria,
% c2_rS1 (type 2 criterion for "S1" responses),
% c2_rS2 (type 2 criterion for "S2" responses)

function [f2 h2] = sdt_2_ot2_roc(c2_rS1, c2_rS2, d, c, pS2)

pS1 = 1 - pS2;

S2mu =  d/2;
S1mu = -d/2;

pHit  = 1 - normcdf(c,S2mu,1);
pCR   = normcdf(c,S1mu,1);

pC = pS2*pHit + pS1*pCR;
pI = 1 - pC;

% "S1"
cdfS1 = normcdf(c2_rS1,S1mu,1);
cdfS2 = normcdf(c2_rS1,S2mu,1);

h2_1 = cdfS1 .* pS1 ./ pC;
f2_1 = cdfS2 .* pS2 ./ pI;

% "S2"
cdfS1 = 1 - normcdf(c2_rS2,S1mu,1);
cdfS2 = 1 - normcdf(c2_rS2,S2mu,1);

h2_2 = cdfS2 .* pS2 ./ pC;
f2_2 = cdfS1 .* pS1 ./ pI;

h2 = h2_1 + h2_2;
f2 = f2_1 + f2_2;

end