function w = tricube_weight(x,x0,h)
%TRICUBE_WEIGHT Summary of this function goes here
%   Detailed explanation goes here
	w = ((x-x0)./h);
	lgc = abs(w) < 1;
	w(~lgc) = 0;
	w(lgc) = (1 - abs(w(lgc)).^3).^3;
end
