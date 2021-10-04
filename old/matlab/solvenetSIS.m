function solvenetSIS(B,G,D,IDX,nodes,maxIter,deltaT,iniz,pert)

%%%% PARAMETERS %%%%%
%
% solvenetSIS(4.5,5,3,10,100,5000,0.005,5000,.1)
% IDX is the index of the metanode which goes beyond the threshold
% B : beta (transmission rate)
% G : gamma (recovery rate)
% D: diffusion rate
% nodes: number of nodes
% 


clc;
close all;
z=1;
k=zeros(1,nodes);


phiStar = 0*ones(1,nodes);

% Generate the adjacency list of the network
A=PAm(nodes,5);
figure(100)
imagesc(A)
%A=random_graph(nodes,0.2);
for i=1:nodes
    k(i)=sum(A(i,:));
end
L=A-diag(k);
[V,D_aut]=eig(L);
diag(D_aut)

phi=zeros(1,nodes);

% Initial values 
perturb=pert*(rand(1,nodes));
phiOld=phiStar+perturb;
 

% File to store the solution 
fileName=sprintf('solMF_nodes_%3.2e_steps_%3.2e',nodes,maxIter);   
fid = fopen(fileName,'w');  

% Save initial values 
fprintf(fid,'%3.5f ',0.0);
for i=1:nodes
   fprintf(fid,'%3.5f ',phiOld(i));
end
fprintf(fid,'\n');


% Different <k> per metanode

B=B*ones(1,nodes);
B(1,IDX)=-(-G+D_aut(IDX,IDX))+5;
B(1,IDX)
IDX
% IDX1=IDX-1
% IDX2=IDX+1
% B(1,IDX1)=-(-G+D_aut(IDX1,IDX1))+5;
% B(1,IDX2)=-(-G+D_aut(IDX2,IDX2))+5;
B-G*ones(1,nodes)+diag(D_aut)'



% Solve the system 
for t=1:maxIter 
k1=zeros(1,nodes);
k2=zeros(1,nodes);
k3=zeros(1,nodes);
k4=zeros(1,nodes);
    for i=1:nodes
        if(phiOld(i)<0)
            disp('Concentrazione negativa');
        end
    end
    for i=1:nodes
        Somma_phi=L(i,:)*phiOld';
        k1(i)=deltaT*(D*Somma_phi+B(i)*phiOld(i)*(1-phiOld(i))-G*phiOld(i));
        phi_k1(i)=phiOld(i)+k1(i)/2;
    end
    for i=1:nodes
        Somma_phi=L(i,:)*phi_k1';
        k2(i)=deltaT*(D*Somma_phi+B(i)*phi_k1(i)*(1-phi_k1(i))-G*phi_k1(i));
        phi_k2(i)=phiOld(i)+k2(i)/2;
    end
    for i=1:nodes
        Somma_phi=L(i,:)*phi_k2';
        k3(i)=deltaT*(D*Somma_phi+B(i)*phi_k2(i)*(1-phi_k2(i))-G*phi_k2(i));
        phi_k3(i)=phiOld(i)+k3(i);
    end
    for i=1:nodes
        Somma_phi=L(i,:)*phi_k3';
        k4(i)=deltaT*(D*Somma_phi+B(i)*phi_k3(i)*(1-phi_k3(i))-G*phi_k3(i));
        phi(i)=phiOld(i)+(1/6)*(k1(i)+2*k2(i)+2*k3(i)+k4(i));
    end

    if(mod(t,z)==0)
        fprintf(fid,'%3.5f ',t*deltaT);
        for i=1:nodes
            fprintf(fid,'%3.5f ',phi(i));
        end
        fprintf(fid,'\n');
    end
    phiOld=phi;
end
    
fclose(fid);


% Movie
data=load(fileName);
[n,~]=size(data);
n=n-1;
time=data(1:n,1);
data=data(1:n,2:nodes+1);



figure(1)
hold on
title('X_i evolution')
xlabel('Time')
ylabel('Nodes')
plot(time,data)
colorbar
hold off


figure(2)
hold on
title('Linear localization')
xlabel('Nodes')
ylabel('X_i (\infty) vs. V_i (1)')
data(iniz,:)=data(iniz,:)-min(data(iniz,:))*ones(1,nodes);
plot(data(iniz,:)/max(abs(data(iniz,:))),'bo')
plot(abs(real(V(:,IDX)))/max(abs(real(V(:,IDX)))),'r*')
%plot(abs(real(V(:,2)))/max(abs(real(V(:,2)))),'r<')
% plot(abs(real(V(:,3)))/max(abs(real(V(:,3)))),'rs')
%plot(abs(real(V(:,100)))/max(real(V(:,1))),'r+')
hold off

% figure(6)
% imagesc(abs(V))
% colorbar
% colormap('jet')
       
