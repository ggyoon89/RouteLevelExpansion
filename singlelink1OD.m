% Random generation (for comparison)

PHR=struct('p',[0,0,0]);
for a=1:200
    tau=randi([13,50]);
    delta=0.25;
    alpha=1+sqrt(log(2/delta)/2);
    Theta=[3.5,-0.3,0.15];
    Y=zeros(tau,1);
    Pick=zeros(tau,1);
    for j=1:100
        A=[];
        X=[];
        Mu=[];
        Z=[];
        for i=1:tau
            A(i,1)=randi([1,50]);
            X(i,:)=T(A(i,1),:);
            Mu(i,1)=-Theta*X(i,:)';
            Z(i,1)=C(A(i,1),1);
            for k=1:Z(i,1)
                if rand<1/(1+exp(Mu(i,1)))
                    Y(i,1)=1;
                end
            end
        end
        W=[A,X(:,2),Z,Mu];

        W=sortrows(unique(W,'rows'),4,'ascend');
        for i=1:size(W,1)
            Pick(i,4*j-3)=W(i,1);
            Pick(i,4*j-2)=W(i,2);
            Pick(i,4*j-1)=W(i,3);
        end
        passhr(j,1)=0;
        total=0;
        for i=1:size(W,1)
            Pick(i,4*j)=0;
            for k=1:W(i,3)
                if rand<1/(1+exp(W(i,4)))
                    Pick(i,4*j)=Pick(i,4*j)+1;
                end
            end
            total=sum(Pick(1:i,4*j-1));
            if total<100
                passhr(j,1)=passhr(j,1)+Pick(i,4*j-1)*Pick(i,4*j-2);
            else
                passhr(j,1)=passhr(j,1)+(Pick(i,4*j-1)-(total-100))*Pick(i,4*j-2);
                passhr(j,2)=i;
                break
            end
            if total>100                
                passhr(j,2)=i;
                passhr(j,3)=tau;
            else
                passhr(j,2)=0;
                passhr(j,3)=tau;
            end
        end
    end
    PHR(1).p=[PHR(1).p;passhr];
end


%Learning

PHR2=struct();
for j=1:15000
        A=[];
        X=[];
        Z=[];
        tau=randi([13,50]);
        delta=0.25;
        alpha=1+sqrt(log(2/delta)/2);
        Theta=[3.5,-0.3,0.15];
        V=zeros(3,3);
        Y=zeros(tau,1);
        for i=1:tau
            A(i,1)=randi([1,50]);
            X(i,:)=T(A(i,1),:);
            Z(i,1)=C(A(i,1),1);
            Y(i,1)=0;
            for k=1:Z(i,1)
                if rand<1/(1+exp(-Theta*X(i,:)'))
                    Y(i,1)=Y(i,1)+1;
                end
            end
            V=V+X(i,:)'*X(i,:)*Z(i,1);    
        end

        for t=tau+1:2*tau
            F=@(thetat)sngllink1ODfcn(thetat,X,Y,Z,t);
            thetat=fsolve(F,zeros(1,3));
            theta(t,:)=thetat;

            armst=zeros(50,1);
            maxa=-inf;
            for i=1:50
                Xta(:,i)=T(i,:)';
                armst(i,1)=thetat*Xta(:,i)+alpha*sqrt(Xta(:,i)'*V^-1*Xta(:,i));
                if armst(i,1)>maxa
                    maxa=armst(i,1);
                    X(t,:)=Xta(:,i)';
                    arms(t,1)=i;
                    arms(t,2)=C(i,1);
                end
            end
            Y(t,1)=0;
            Z(t,1)=arms(t,2);
            for k=1:Z(t,1)
                if rand<1/(1+exp(-Theta*X(t,:)'))
                    Y(t,1)=Y(t,1)+1;
                end
            end
            V=V+X(t,:)'*X(t,:)*Z(t,1);

        end
        for k=1:50
            armexit(k,1)=k;
        end
        armexit(:,2)=armst;
        armexit(:,3)=C;
        armexit(:,4)=T(:,2);
        armexit(:,5)=T(:,3);
        armexit=sortrows(armexit,2,'descend');
        k=1;
        sm=0;

        while k>0
            sm=sm+armexit(k,3);
            if sm>100
                chosen(1:k,1:5)=armexit(1:k,1:5);            
                break
            end
            k=k+1;
        end
    PHR2.passhr(j,3)=tau;
    PHR2.passhr(j,2)=k;
    PHR2.passhr(j,1)=sum(chosen(1:k,4).*chosen(1:k,3))-chosen(k,4)*(sm-100);
    PHR2.thetat(j,:)=thetat;
end