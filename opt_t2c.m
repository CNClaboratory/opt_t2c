function opt = opt_t2c(type, varargin)
% opt = opt_t2c(type, varargin)
%
% Compute the type 2 criteria that optimize a given target measure according
% to classical SDT, given relevant factors such as d', c, p(S2), etc.
%
% TYPE
% ----
% The input 'type' is a string that determines what target measure is to be
% optimized. Possible values of type are:
%
% * 't2acc'       - p(type 2 correct)
% * 't2reward'    - expected reward based on type 2 reward contingency table
% * 'calibration' - reporting high confidence for trials where p(correct) >
%                   some p(correct) threshold
% * 'hr2-far2'    - difference b/t type 2 hit rate and type 2 false alarm rate
%
% INPUTS BY TYPE
% --------------
% Use the following input styles for each type:
%
% opt = opt_t2c('t2acc', d, c, pS2)
% opt = opt_t2c('t2reward', d, c, pS2, Q2)
% opt = opt_t2c('calibration', d, c, pS2, OT)
% opt = opt_t2c('hr2-far2', d, c)
%
% where 
% * d   - SDT measure of sensitivity d'
% * c   - SDT measure of response bias c
% * pS2 - probability of S2 stimulus being presented
% * Q2  - the relative reward quotient (R_CR2 - R_FA2) / (R_hit2 - R_miss2),
%         where R values denote the amount of points earned following each 
%         type 2 event
% * OT  - the odds ratio p(correct)_T / (1-p(correct)_T), where p(correct)_T 
%         is the threshold level of p(correct) such that confidence is high 
%         when estimated p(correct) > p(correct)_T
%
% OUTPUT
% ------
% The output struct "opt" contains the following fields:
% * opt.c2_rS1 - optimal type 2 criterion for "S1" responses
% * opt.c2_rS2 - optimal type 2 criterion for "S2" responses
% * opt.B2_rS1 - optimal type 2 beta for "S1" responses
% * opt.B2_rS1 - optimal type 2 beta for "S2" responses
% * opt.input  - a copy of the input variables used to create "opt"
%
% 3/19/2022 Brian Maniscalco, Lucie Charles, & Megan Peters

switch type
    
    % optimizing type 2 accuracy    
    case 't2acc'
        
        % parse input
        d   = varargin{1};
        c   = varargin{2};
        pS2 = varargin{3};
        
        pS1 = 1 - pS2;
        B   = exp(c .* d);

        % compute optimal criteria
        opt.c2_rS1 = min( log(pS1./pS2) ./ d, c );
        opt.c2_rS2 = max( log(pS1./pS2) ./ d, c );
        opt.B2_rS1 = min( pS1./pS2, B );        
        opt.B2_rS2 = max( pS1./pS2, B );

        % package input info into output
        opt.input.type = type;
        opt.input.d    = d;
        opt.input.c    = c;
        opt.input.pS2  = pS2;
        
    % optimizing type 2 reward
    case 't2reward'

        % parse input
        d   = varargin{1};
        c   = varargin{2};
        pS2 = varargin{3};
        Q2  = varargin{4};
        
        pS1 = 1 - pS2;
        B   = exp(c .* d);

        % compute optimal criteria
        opt.c2_rS1 = min( ( log(pS1./pS2) - log(Q2) ) ./ d, c);
        opt.c2_rS2 = max( ( log(pS1./pS2) + log(Q2) ) ./ d, c);
        opt.B2_rS1 = min( (pS1./pS2) .* Q2, B);
        opt.B2_rS2 = max( (pS1./pS2) .* Q2, B);

        % package input info into output
        opt.input.type = type;
        opt.input.d    = d;
        opt.input.c    = c;
        opt.input.pS2  = pS2;
        opt.input.Q2   = Q2;
        
    % calibrating confidence to expected accuracy
    case 'calibration'
        
        % parse input
        d   = varargin{1};
        c   = varargin{2};
        pS2 = varargin{3};
        OT  = varargin{4};
        
        pS1 = 1 - pS2;
        B   = exp(c .* d);

        % compute optimal criteria
        opt.c2_rS1 = min( ( log(pS1./pS2) - log( OT ) ) ./ d, c);
        opt.c2_rS2 = max( ( log(pS1./pS2) + log( OT ) ) ./ d, c);
        opt.B2_rS1 = min( (pS1./pS2) .* (1/OT), B);
        opt.B2_rS2 = max( (pS1./pS2) .* OT, B);

        % package input info into output
        opt.input.type = type;
        opt.input.d    = d;
        opt.input.c    = c;
        opt.input.pS2  = pS2;
        opt.input.OT   = OT;
        
    % type 2 HR - FAR
    case 'hr2-far2'
        
        % parse input
        d   = varargin{1};
        c   = varargin{2};
        
        HR  = 1 - normcdf(c, d/2);  % hit rate
        FAR = 1 - normcdf(c, -d/2); % false alarm rate
        B   = exp(c .* d);

        % compute optimal criteria
        opt.c2_rS1 = min( log( (1-HR) ./ (1-FAR) ) ./ d, c );
        opt.c2_rS2 = max( log( HR./FAR ) ./ d, c );
        opt.B2_rS1 = min( (1-HR) ./ (1-FAR), B );
        opt.B2_rS2 = max( HR ./ FAR, B );

        % package input info into output
        opt.input.type = type;
        opt.input.d    = d;
        opt.input.c    = c;
        
    otherwise
        error('invalid type entered')

end