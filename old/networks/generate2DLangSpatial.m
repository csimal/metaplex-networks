%function A = generate1DSpatialNet(numNodes, alpha, lambda)
clear all;
close all;
%% Change these to tune the network
numNodes = 1000;
k = 12;
sigma = 0.3;
expectedDegree = poissrnd(k,numNodes,1);

%% Generate the network
coords = rand(numNodes,2);
eulerDist = dist(coords');
k = sum(expectedDegree)/numNodes;
density = numNodes;

% Distance kernel
f = exp((-eulerDist.^2)/(2*sigma^2))/(2*pi*sigma^2);

linkProbabilites = zeros(numNodes,numNodes);
for i = 1:numNodes
    for j = i:numNodes
        linkProbabilites(i,j) = expectedDegree(i)*expectedDegree(j)*f(i,j)/(density*k);
    end
end

% Randomly keep edges
A = linkProbabilites>rand(numNodes,numNodes);

%Remove self connections and ensure symmetric
A = A-diag(diag(A));
A = A + A';
%A = A>0;

[M,I] = find(sum(A,2)==0);
M = flipud(M);
for i = M
    A(i,:) = [];
    A(:,i) = [];
    coords(i,:) = [];
end

%% Plotting
imagesc(A)

g = graph(A)
figure
h1 = plot(g);
h1.XData = coords(:,1);
h1.YData = coords(:,2);

figure
h2 = plot(g);
h2.XData = coords(:,1);
h2.YData = coords(:,2);
nodeIDs1 = nearest(g,1,1);
nodeIDs2 = nearest(g,2,1);
nodeIDs3 = nearest(g,3,1);
highlight(h2,[1;nodeIDs1],'Nodecolor','g');
highlight(h2,[2;nodeIDs2],'Nodecolor','r');
highlight(h2,[3;nodeIDs3],'Nodecolor','m');
h2.MarkerSize = 2;
h2.EdgeColor = [1,1,1];
highlight(h2,[1,2,3],'MarkerSize',10);
highlight(h2,[nodeIDs1;nodeIDs2;nodeIDs3],'MarkerSize',6);


isConnected(A)
%end