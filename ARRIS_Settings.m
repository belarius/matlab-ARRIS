function bp = ARRIS_Settings(varargin)
%ARRIS_SETTINGS Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(varargin)
    if iscell(varargin{1})
        varargin = varargin{1};
    end
end

bp = struct();
bp.burn_iter = 4000;		%		Burn-in period
bp.samp_iter = 16000;		%		Sample MCMC interations
bp.k = 150;					%		Initial number of knots
bp.beta_iter = 6;			%		Iterations for parameter resampling
bp.probbd = .4;				%		Birth/Death probability
bp.tau = 2;					%		Tension in the beta distribution sampler (higher = stickier prior)
bp.conf_level = .95;		%		Alpha level
bp.threshold = -5;			%		Late-stage rejection likelihood threshold
bp.prior_id = 'UNIFORM';	%		Prior distribution of knots
bp.dparams = 6;				%		Non-uniform prior distribution parameter(s)
bp.iparams = [1 150];		%		Uniform prior distribution parameters
bp.range = [0 5000];		%		Range sampled
bp.bin_width = 10;			%		Window size when data are discretized
bp.max_knots = 200;			%		Upper limit on model complexity


if mod(length(varargin),2)
    error('CPRBayes:BadlyFormedParameters','Incorrectly formatted parameters');
else
    rsargin = reshape(varargin,2,length(varargin)./2)';
    for i = 1:length(rsargin(:,1))
        switch rsargin{i,1}
            case 'burn_iterations'
                bp.burn_iter = rsargin{i,2};
            case 'samp_iterations'
                bp.samp_iter = rsargin{i,2};
            case 'init_knots'
                bp.k = rsargin{i,2};
            case 'beta_iterations'
                bp.beta_iter = rsargin{i,2};
            case 'prob_birth_death'
                bp.probbd = rsargin{i,2};
            case 'tau'
                bp.tau = rsargin{i,2};
            case 'conf_level'
                bp.conf_level = rsargin{i,2};
            case 'threshold'
                bp.threshold = rsargin{i,2};
            case 'prior_id'
                bp.prior_id = rsargin{i,2};
            case 'dparams'
                bp.dparams = rsargin{i,2};
            case 'iparams'
                bp.iparams = rsargin{i,2};
            case 'range'
                bp.range = rsargin{i,2};
            case 'bin_width'
                bp.bin_width = rsargin{i,2};
            case 'max_knots'
                bp.max_knots = rsargin{i,2};
        end
    end
end


end

