function [GCV,qmap,dmap,bw,spk_struct] = NeuroKD_GCV_single(D,t,t_bins,h,r_index,bw,spk_struct)
%NEUROKD_GCV_SINGLE Summary of this function goes here
%   Detailed explanation goes here

tlen = length(t_bins);
dlen = length(D);

if nargin < 7
	spk_struct = struct('spk',[],'spk_d',[],'spk_l',[]);
	spk_struct.spk = cell2mat(D);
	spk_struct.spk_d = D;
	spk_struct.spk_l = zeros(dlen,1);
	for i = 1:dlen
		spk_struct.spk_d{i}(:) = i;
		spk_struct.spk_l(i) = length(D{i});
	end
	spk_struct.spk_d = cell2mat(spk_struct.spk_d);
end

if nargin < 6
	bw = zeros(dlen,1);
	for i = 1:dlen
		bw(i) = 1.06.*((median(diff(D{i})))./log(2)).^(6/5);
	end
end

if nargin < 5
	r_index = (1:dlen)';
end

if h > 1
	d = D{t};
	kd = ksdensity(d,t_bins,'bandwidth',bw(t),'kernel','box');
	dmap = kd.*length(d)./sum(kd);

	d = [spk_struct.spk zeros(size(spk_struct.spk,1),1)];
	tcw = tricube_weight(r_index,r_index(t),h);
	w = sum(tcw);
	d(:,2) = tcw(spk_struct.spk_d);
	dex = (d(:,2)>0);
	
	if sum(dex)==0
		qmap = NaN;
		GCV = NaN;
	else
%		kd = local_box_ksd(d(dex,1),t_bins,tlen,bw,d(dex,2));
		kd = ksdensity(d(dex,1),t_bins,'bandwidth',bw(t),'weights',d(dex,2),'kernel','box');
		qmap = (kd.*sum(d(dex,2))./(w.*sum(kd)))';
		N = tlen;
		Trc = tlen./w;
		SS = sum((dmap(:)-qmap(:)).^2);
		GCV = N*SS/((N-1.5*Trc)^2);
	end
else
	d = D{t};
%	kd = local_box_ksd(d,t_bins,tlen,bw,ones(length(d),1));
	kd = ksdensity(d,t_bins,'bandwidth',bw(t),'kernel','box');
	dmap = kd.*length(d)./sum(kd);
	qmap = dmap;
	GCV = NaN;
end

end

function kd = local_box_ksd(d,t_bins,tlen,bw,w)
	bw = bw.*3.464101615137754; %scaled by sqrt(12)
	M = abs(repmat(d,1,tlen)-repmat(t_bins,length(d),1));
	M = (M<(bw./2))./(bw);
	kd = M'*w;
end
