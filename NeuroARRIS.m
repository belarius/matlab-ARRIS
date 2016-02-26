function [fit, stats] = NeuroARRIS(data,trials,varargin)
%NEUROARRIS Adaptive Rate Regression Involving Steps for
%Poisson-distributed spike rasters.
%   Detailed explanation goes here

%====================
%===INITIALIZATION===
%====================

fit = 0;
stats = NaN;

%==Parameters==
if nargin < 2
    trials = 1;
end

if nargin < 3
	varargin = [];
end
bp = ARRIS_Settings(varargin);

%==Data Setup==
if isempty(data)
	return;
else
	if size(data,2) == 1
		data = [data ones(length(data),1)];
	elseif size(data,2) > 2
		%ERROR MESSAGE
	end
	bd = struct();
	bd.range = bp.range;
	bd.x_bins = bd.range(1):bp.bin_width:bd.range(2);
	bd.y = zeros(length(bd.x_bins),1);
	bd.n = length(bd.x_bins);
	for i = 1:size(data,1)
		dex = find(bd.x_bins==floor(data(i,1)/bp.bin_width)*bp.bin_width,1);
		bd.y(dex) = bd.y(dex) + data(i,2);
	end
	bd.xg = (bd.x_bins - bd.range(1) + bp.bin_width/2)/(bd.range(2) - bd.range(1) + bp.bin_width);
	bd.scale_factor = 1./(trials);
	bd.trials = trials;
end

%==Model Setup==
m = struct();
m.n = bd.n;
m.fit_info = 0;
m.accept_prob = 0;
m.knots_interior = (1:bp.k)'/(bp.k+1);

%=====================
%====STARTER MODEL====
%=====================

%==Mu Setup==
mu_start = bd.y;
bp.twist = sum(mu_start);
mu_start(mu_start == 0) = 0.5/length(mu_start);
bp.twist = sum(mu_start)./bp.twist;

%==Initial Fit Attempt==
m = setBasisAndFitModel(m,bd,mu_start);

%==MCMC Setup==
[birth_probs,death_probs] = getRJProbs(bp);		% using the user parameters, build the knot prior
tot_iter = bp.samp_iter + bp.burn_iter;			% total iterations is the sample iterations + the burn in iterations
models = cell(bp.samp_iter,1);					% models is the output - the collection of fitted models after burn-in
moves = zeros(3,1);

%==MCMC==
i = 1;
j = 0;
while i <= tot_iter
	j = j+1;
	u = rand(1);
	k = length(m.knots_interior);
	if u < birth_probs(k)
		bdr = 1; % birth step
		mNew = getBirthModel(m,bp,bd,mu_start);
	elseif 1-u < death_probs(k)
		bdr = 2; % death step
		mNew = getDeathModel(m,bp,bd,mu_start);
	else
		bdr = 3; % relocation step
		mNew = getRelocateModel(m,bp,bd,mu_start);
	end
	% Check is two knots are too close together to include any data
	dex = find(sum(mNew.basis(:,3:end))==0);
	while ~isempty(dex)
		if dex(end) > length(mNew.knots_interior)
			dex(end) = dex(end) - 1;
		elseif dex(end) > 1
			dex(end) = randsample([dex(end)-1;dex(end)],1);
		end
		mNew.knots_interior(dex(end)) = [];
		mNew.knots_all(dex(end)+4) = [];
		mNew.basis(:,dex(end)+2) = [];
		dex(end) = [];
		mNew = getBirthModel(mNew,bp,bd,mu_start);
		dex = find(sum(mNew.basis(:,3:end))==0);
	end
	moves(bdr) = moves(bdr) + 1;
	if rand(1) < mNew.accept_prob
		m = mNew;
	end
	if i > bp.burn_iter
		[m,laterej] = RandMu(m,bd,bp);
		if laterej
			i = i - 1;
		else
			models{i-bp.burn_iter} = getStats(m,bd);
		end
	end
	i = i + 1;
end

%==Output Prep==
models = cell2mat(models);
modelCell = struct2cell(models);

%==Output==
fit = mean(cell2mat(modelCell(4,:)).*bd.scale_factor./bp.twist,2);
stats = cell2mat(modelCell(4,:)).*bd.scale_factor./bp.twist;

end

%=================
%====FUNCTIONS====
%=================

function [m,rcond] = setBasisAndFitModel(m,bd,mu_start)  
% Description
	%Pad extra knots onto the ends
	m.knots_all = [zeros(4,1); m.knots_interior; ones(4,1)];
	m.basis = getModelSegments(m, bd);
    [m,rcond] = FitPoissonModel(m,bd,mu_start);
end

function basis = getModelSegments(m, bd)
% return a new model basis for the knots in m

	%Set variables
	knots = m.knots_all;
	knots_int = m.knots_interior;
	l_knots_int = length(knots_int)+3;
	xg = bd.xg;
	l_xg = length(xg);
	basis = ones(l_xg,l_knots_int);
	for i = 3:l_knots_int
		basis(:,i) = (xg' > knots(i+1)) .* (xg' <= knots(i+2));
	end
end

function [m,rcond] = FitPoissonModel(m,bd,mu_start)
% given a set of knots and a precomputed basis, fit a poisson model
	warning('off','MATLAB:nearlySingularMatrix')
	warning('off','MATLAB:singularMatrix')
	p = length(m.knots_interior);
	x = m.basis(:,3:end);
	y = bd.y;
	B = zeros(p+1,1);
	for i = 1:p+1;
		B(i) = sum(mu_start.*x(:,i))./sum(x(:,i));
		if B(i) == 0
			B(i) = 0.5./length(y);
		end
		if isnan(B(i))
			B(i) = 0.5./length(y);
		end
	end
	mu = x*B;
 	loglik1 = sum(y .* log(mu) - mu);
	h = diag(mu_start) * x;
	J = h'*x;
	[J_chol,p] = chol(J);
 	if p > 0
		converged = 0;
		rcond = -1;
		J = diag(B)*(trace(J)./sum(B));
		[J_chol,p] = chol(J);
		J = [J_chol B];
	else
		converged = 1;
		rcond = NaN;
		J = [J_chol B];
	end
	m.loglikelihood = loglik1;
	m.sigma = 1;
	m.fit_info = converged;
	m.J = J;
	warning('on','MATLAB:nearlySingularMatrix')
	warning('on','MATLAB:singularMatrix')
end

function [birth,death] = getRJProbs(bp)
% builds the knot prior
	i = 1:(bp.max_knots-1);
	if strcmp(bp.prior_id,'POISSON')
		pratio = bp.dparams(1) ./ (i+1);
	elseif strcmp(bp.prior_id,'UNIFORM')
		pratio = (i >= bp.iparams(1)) & (i < bp.iparams(2));
	elseif strcmp(bp.prior_id,'USER')
		pratio = i*0;
		pratio(bp.iparams(1):bp.iparams(2)-1) = bp.dparams(i+1)./bp.dparams(i);
	else
		pratio = i*0;
	end
	pratio = [0 pratio 0];
	birth = zeros(bp.max_knots+1,1);
	death = zeros(bp.max_knots+1,1);
	c = bp.probbd;
	birth(1) = 0;
	death(1) = 0;
	birth(2) = (pratio(2) >= eps) * c * min(1,pratio(2));
	death(2) = 0;
	if pratio(end-1) < eps
		death(end) = c;
	else
		death(end) = c * min(1,1/pratio(end-1));
	end
	for i=3:bp.max_knots
		birth(i) = (pratio(i) >= eps) * c * min(1,pratio(i));
		if pratio(i-1) < eps
			death(i) = c;
		else
			death(i) = c * min(1,1/pratio(i-1));
		end
	end
end

function m = getBirthModel(m,bp,bd,mu_start)
% compute a new model, adding a single knot
	oldLike = m.loglikelihood;
	k = length(m.knots_interior);
	cand = getBirthKnotInfo(m.knots_interior,bp.tau);
	dens = getTransitionDensity(cand,bp.tau,m.knots_interior);
	m.knots_interior= sort([m.knots_interior; cand]);
	m = setBasisAndFitModel(m,bd,mu_start);
	if m.fit_info > 0
		m.accept_prob = exp(m.loglikelihood - oldLike + log(k) - log(dens) - .5 * log(m.n));
	else
		m.accept_prob = 0;
	end
end

function m = getDeathModel(m,bp,bd,mu_start)
% compute a new model, deleting a knot
	oldLike = m.loglikelihood;
	k = length(m.knots_interior);
	dkPos = randi(k);
	cand = m.knots_interior(dkPos);
	m.knots_interior(dkPos) = [];
	dens = getTransitionDensity(cand,bp.tau,m.knots_interior);
	m = setBasisAndFitModel(m,bd,mu_start);
	if m.fit_info > 0
		m.accept_prob = exp(m.loglikelihood - oldLike - log(k-1) + log(dens) + .5 * log(m.n));
	else
		m.accept_prob = 0;
	end
end

function m = getRelocateModel(m,bp,bd,mu_start)
% compute a new model, relocating a single knot
	oldLike = m.loglikelihood;
	[birth_cand,dkPos] = getBirthKnotInfo(m.knots_interior,bp.tau);
	death_cand = m.knots_interior(dkPos);
	m.knots_interior(dkPos) = [];
	m.knots_interior = sort([m.knots_interior; birth_cand]);
	dens1 = dbeta(birth_cand,bp.tau*death_cand,bp.tau*(1-death_cand));
	dens2 = dbeta(death_cand,bp.tau*birth_cand,bp.tau*(1-birth_cand));
	m = setBasisAndFitModel(m,bd,mu_start);
	if m.fit_info > 0
		m.accept_prob = exp(m.loglikelihood - oldLike + log(dens2) - log(dens1));
	else
		m.accept_prob = 0;
	end
end
       
function d = getTransitionDensity(cand,tau,knots)
	d = 0;
	for i = 1:length(knots)
		d = d + dbeta(cand,tau*knots(i),tau*(1-knots(i)));
	end
end

function d = dbeta(x,a,b)
	if x <= 0 || x >= 1
		d = 0;
	else
		d = exp((a-1)*log(x) + (b-1)*log(1-x) - betaln(a,b));
	end
end

function [cand,k_rand] = getBirthKnotInfo(knots,tau)
% figure out where a new knot would be added
	k = length(knots);
	k_rand = randi(k,1);
	alpha = tau * knots(k_rand);
	beta = tau - alpha;
	cand = betarnd(alpha,beta);
	while check_cand(cand,knots)
		cand = betarnd(alpha,beta);
	end
end

function iout = check_cand(cand,knots)
	iout = 0;
	iout = iout + (cand < eps);
	iout = iout + (cand > 1-eps);
	for i=1:length(knots)
		iout = iout + (abs(cand-knots(i)) < eps);
	end
end

function [m,laterej] = RandMu(m,bd,bp)
% compute a sample mu, and also return whether or not there should be a late rejection
	[beta,laterej] = RandBeta(m,bd,bp);
	m.random_mu = exp(m.basis(:,3:end) * beta);
end

function [beta,laterej] = RandBeta(m,bd,bp)
% compute a perturbed beta vector
	MHiters = bp.beta_iter;
	threshold = bp.threshold;
	p = length(m.knots_interior)+1;
	lastbeta = m.J(:,end);
	m.random_mu = m.basis(:,3:end) * lastbeta;
	lastfi = sum(bd.y .* log(m.random_mu) - m.random_mu);
	lastgi = sum((m.J(2:end,2:end-1)*lastbeta(2:end)).^2)/(-2*m.n);
	lasthi = 0;
	countMH = 0;
	MHi = 1;
	while MHi <= MHiters
		curbeta = randn(p,1);
		curhi = sum(curbeta.^2);
		curbeta = m.J(:,1:end-1)\curbeta;
		curhi = curhi * -.5;
		curbeta = curbeta + m.J(:,end);
		m.random_mu = m.basis(:,3:end) * curbeta;
		curfi = bd.y .* log(m.random_mu) - m.random_mu;
		curfi = sum(curfi(~isnan(curfi)));
		curgi = sum((m.J(2:end,2:end-1)*curbeta(2:end)).^2)/(-2*m.n);
		r = (curfi - lastfi) + (curgi - lastgi) - (curhi - lasthi);
		if r > 0
			r = 0;
		end
		u = rand(1);
		if u < eps
			u = r-1;
		else
			u = log(u);
		end
		if MHi == 1 && r > threshold
			MHi = MHiters+1;
			u = r - 1;
		end
		if u < r
			lastbeta = curbeta;
			lastfi = curfi;
			lastgi = curgi;
			lasthi = curhi;
			countMH = countMH + 1;
		end	
		MHi = MHi + 1;
	end
	beta = lastbeta;
	laterej = 0;%countMH == 0;
end

function output = getStats(m,bd)
% compute some final statistics for this model m
	output = struct();
	[maxY,maxInd] = max(m.random_mu);
	k = length(m.knots_interior);
	w = k + 1.5;
	v = sum(bd.y .* log(m.random_mu) - m.random_mu);
	output.BIC = v - w*log(m.n);
	output.logLikelihood = v;
	output.knots = m.knots_interior;
	output.sampMu = m.basis(:,3:end) * m.J(:,end);%m.random_mu;
end


