%function A = generate1DSpatialNet(numNodes, alpha, lambda)
numNodes = 484;
n = sqrt(numNodes);
alpha =0.2;
lambda = 0.4;
coords = [];
for i = 1:n
    for j = 1:n
        coords = [coords;[i,j]];
    end
end
        

eulerDist = dist(coords');

linkProbabilites = zeros(numNodes,numNodes);
linkProbabilites = alpha*lambda.^(eulerDist-1);

A = linkProbabilites>rand(numNodes,numNodes);

A = A-diag(diag(A));
A = A + A';
A = A>0;
imagesc(A)

g = graph(A)
figure
h1 = plot(g);
h1.XData = coords(:,1);
h1.YData = coords(:,2);
figure
h2 = plot(g);

isConnected(A)
%end