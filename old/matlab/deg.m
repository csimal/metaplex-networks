function [dgi,dgo]=deg(Adj)
% function [dgi,dgo]=deg(Adj)
%
% This function computes the (non normalized) ingoing and outgoing degree of a network 
% defined via its adjacency matrix.
%
% Timoteo Carletti 23/09/2008
%
% input : Adj Adjacency matrix
%
% output : dgi ingoing degree
%        : dgo outgoing degree

dgi=degin(Adj);
dgo=degout(Adj);