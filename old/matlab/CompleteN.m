function Adj=CompleteN(N)
% function Adj=CompleteN(N)
%
% This file creates the complete networks with N nodes
%
% Timoteo Carletti 7/11/2014
%
%
% input : N number of nodes
%
% output : Adj Adjacency matrix
%

  v=ones(1,N);
  Adj=ones(N,N)-diag(v);

end
