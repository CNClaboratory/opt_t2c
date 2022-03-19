function x2 = type2noisySignalLoss(x1, param)
% x2 = type2noisySignalLoss(x1, param)
%
% Given a 1xN vector of type 1 evidence samples x1, returns a 1xN vector of
% type 2 evidence samples x2 where x2 is subject to noise and signal loss
% according to the parameters param.sd2 and param.k, respectively. More
% specifically,
%
% x2 = (1-param.k)*x1 + normrnd(0, param.sd2, size(x1));

x2 = (1-param.k)*x1 + normrnd(0, param.sd2, size(x1));
