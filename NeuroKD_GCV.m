function [GCV,qmap,dmap,bw,spk_struct] = NeuroKD_GCV(D,t_bins,h,dmap,r_index,bw,spk_struct)
%NEUROKD_GCV Summary of this function goes here
%   Detailed explanation goes here

tlen = length(t_bins);
dlen = length(D);
GCV = 0;

N = dlen.*tlen;
Trc = 0;

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
		b = gamfit(diff(D{i}));
		bm = b(1).*b(2);
		bsd = sqrt(b(1).*b(2).^2);
		bw(i) = 1.06*bsd*(bm)^0.2;
	end
end

if nargin < 5
	r_index = (1:dlen)';
end

if nargin < 4
	dmap = [];
end

if isempty(dmap)
	dmap = zeros(dlen,tlen);
	for j = 1:dlen
%		kd = local_box_ksd(D{j},t_bins,tlen,bw(j),ones(length(D{j}),1));
		kd = ksdensity(D{j},t_bins,'bandwidth',bw(j),'kernel','box');
		dmap(j,:) = kd.*spk_struct.spk_l(j)./sum(kd);
	end
end

if length(h) == 1
	if h > 1
		h = repmat(h,dlen,1);
	else
		qmap = dmap;
		GCV = NaN;
	end
end

if ~isnan(GCV)
	qmap = zeros(dlen,tlen);
	d = [spk_struct.spk zeros(size(spk_struct.spk,1),1)];
	for j = 1:dlen
		tcw = tricube_weight(r_index,r_index(j),h(j));
		w = sum(tcw);
		d(:,2) = tcw(spk_struct.spk_d);
		dex = (d(:,2)>0);
%		kd = local_box_ksd(d(dex,1),t_bins,tlen,bw(j),d(dex,2));
		kd = ksdensity(d(dex,1),t_bins,'bandwidth',bw(j),'weights',d(dex,2),'kernel','box');
		qmap(j,:) = (kd.*sum(d(dex,2))./(w*sum(kd)))';
		Trc = Trc + tlen/w;
	end
	SS = sum((dmap(:)-qmap(:)).^2);
	GCV = N*SS/((N-1.5*Trc)^2);
end
trash=1;
end

function kd = local_box_ksd(d,t_bins,tlen,bw,w)
	bw = bw*3.464101615137754; %scaled by sqrt(12)
	M = abs(repmat(d,1,tlen)-repmat(t_bins,length(d),1));
	M = (M<(bw/2))./(bw);
	kd = M'*w;
end
