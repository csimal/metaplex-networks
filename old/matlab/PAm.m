function Adj=PAm(N,m)
% function Adj=PAm(N,m)
%
% This function creates a undirected Scale Free Network with N vertices according to the standard
% preferential attachement scheme (i.e. gamma=3) adding each time step 1
% new node connected with m already existing nodes proportionally with their degree
%
% Timoteo Carletti 03/01/2015
%
% input : N number of nodes
%         m : number of new links done each time step
%
% output : Adj Adjacency matrix
%
% required functions:
% CompleteN.m
% deg.m

  Adj=zeros(N);

  % the nodes 1...m form a complete graph
  Adj(1:m,1:m)=CompleteN(m);
  
  tt=m+1;
  while (tt<=N)
    [dgi,~]=deg(Adj);
    Pi=cumsum(dgi(1:tt-1))./sum(dgi(1:tt-1));
    
    dL=0;
    L=[];
    while (dL<m)
        r=rand(1);
        kk=length(find(Pi<r))+1;
        L=union(L,kk);
        dL=length(L);
    end
    % add a new node with m links to the already existing nodes
    Adj(tt,L)=ones(1,dL);
    Adj(L,tt)=ones(dL,1);
    %
    tt=tt+1;  
  end

end

