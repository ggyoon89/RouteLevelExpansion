function [Obj] = totdem(ExpThetaX,NodeInfo,Wk)
SumExp=Wk*ExpThetaX;

for i=1:size(ExpThetaX,2)
    B(1,i)=NodeInfo(i,2)*SumExp(1,i)/(1+SumExp(1,i));
end
    Obj=-sum(B); % minimize negative sum
end