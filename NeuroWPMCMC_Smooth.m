function [activityMap,activityVar,activityDst] = NeuroWPMCMC_Smooth(D,t_bins,r_fit,h,prctls,r_index,dist_approx_range)
%NEUROWPMCMC_SMOOTH Summary of this function goes here
%   Detailed explanation goes here

if nargin < 7
	dist_approx_range = 100;
end

if nargin < 6
	r_index = (1:length(D))';
end

if nargin < 5
	prctls = [];
end

dlen = length(r_index);
tlen = length(t_bins);
rlen = length(r_fit);
plen = length(prctls);

if length(h) == 1
	h = repmat(h,rlen,1);
end

activityMap = zeros(rlen,tlen);
activityVar = zeros(rlen,tlen,plen);
activityDst = zeros(rlen,tlen,dist_approx_range);

for j = 1:rlen
	d = [];
	w = 0;
	for i = 1:dlen
		d = [d;D{i} ones(length(D{i}),1).*tricube_weight(r_index(i),r_fit(j),h(j))];
		w = w + tricube_weight(r_index(i),r_fit(j),h(j));
	end
	[q,q_d] = NeuroARRIS(d,w);
	activityMap(j,:) = q';
	[qlen,ilen] = size(q_d);
	for i = 1:qlen
		q_d(i,:) = sort(q_d(i,:));
	end
	if plen > 0
		for p = 1:plen
			activityVar(j,:,p) = q_d(:,ilen.*prctls(p)/100);
		end
	end
	if dist_approx_range > 0
		for p = 1:100
			activityDst(j,:,p) = q_d(:,ilen.*(p-0.5)/100);
		end
	end
end


end


