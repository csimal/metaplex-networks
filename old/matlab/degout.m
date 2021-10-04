function degree=degout(Adj)
% function degree=degout(Adj)
%
% This function computes the (non normalized) outgoing degree of a network defined via its adjacency matrix.
%
% Timoteo Carletti 23/09/2008
%
% input : Adj Adjacency matrix
%
% output : degree outgoing degree
%

degree=sum(Adj,2);
