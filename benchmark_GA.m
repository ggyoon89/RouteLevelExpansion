load('casestudy2.mat')
% x = ga(@totdem,nvars,A,b,[],[],lb,ub,nonlcon,IntCon);
%options=optimoptions('ga','ConstraintTolerance',0,'FunctionTolerance',1e-8,'MaxGeneration',60000,'PopulationSize',10000)
G=@(Wk)totdem(ExpThetaX,NodeInfo,Wk);
% [x,fval,exitflag,output,population,scores] = ga(G,60,A,b,[],[],[],[],[],(1:60))

t=1000;
RT=zeros(61,t);
for i=1:t
    x=ga(G,60,A,b,[],[],lb,ub,[],(1:60));
    RT(2:61,i)=x';
    RT(1,i)=-G(x);
    i
end
[M1,M2]=max(RT(1,:)');