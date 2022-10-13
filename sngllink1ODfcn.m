function F=sngllink1ODfcn(thetat,X,Y,Z,t)

F=zeros(1,3);
for i=1:t-1
    F=F+(Y(i,1)-Z(i,1)*(1+exp(-thetat*X(i,:)'))^-1)*X(i,:);
end