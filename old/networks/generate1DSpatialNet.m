%function A = generate1DSpatialNet(numNodes, alpha, lambda)
numNodes = 100;
alpha = 1;
lambda = 0.65;

coords = [1:numNodes];%rand(1,numNodes);
dist = zeros(numNodes,numNodes);
for i = 1:numNodes
    dist(i,:) = sqrt((coords(i)-coords).^2);
end

linkProbabilites = zeros(numNodes,numNodes);
linkProbabilites = alpha*lambda.^(dist-1);

A = linkProbabilites>rand(numNodes,numNodes);

A = A-diag(diag(A));
A = A + A';
A = A>0;
imagesc(A)
%end