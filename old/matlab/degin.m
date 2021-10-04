function degree=degin(Adj)
% function degree=degin(Adj)
%
% This function computes the (non normalized) ingoing degree of a network defined via its adjacency matrix.
%
% Timoteo Carletti 23/09/2008
%
% input : Adj Adjacency matrix
%
% output : degree ingoing degree
%

degree=sum(Adj);
