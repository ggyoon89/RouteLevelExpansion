% % Route enumeration
% % Extract node information from link node incidence
% 
% Conn=struct();
% for a=1:25
%     Conn(a).Next=[];
%     Conn(a).TT=[];
%     Conn(a).Cap=[];
%     Conn(a).Ind=[];
%     for b=1:24
%         if LinkNode(b,2) == a
%             Conn(a).Next=[Conn(a).Next;LinkNode(b,3)];
%             Conn(a).TT=[Conn(a).TT;LinkNode(b,5)];
%             Conn(a).Cap=[Conn(a).Cap;LinkNode(b,4)];
%             Conn(a).Ind=[Conn(a).Ind;LinkNode(b,1)];
%         end
%     end
% end

% % Routes with 1 link
% n=1;
% Enumeroute=struct;
% Nodei=13;
% Btt=0;
% Blink=0;
% K(n)=1;
% for i=1:size(Conn(Nodei).Next,1)
%     Broute=[Nodei,Conn(Nodei).Next(i,1)];
%     Btt=[0,Conn(Nodei).TT(i,1)];
%     Bcap=[0,Conn(Nodei).Cap(i,1)];
%     Blink=[0,Conn(Nodei).Ind(i,1)];
%     Enumeroute(n).Route(K(n),:)=Broute;
%     Enumeroute(n).TT(K(n),:)=Btt;
%     Enumeroute(n).Cap(K(n),:)=Bcap;
%     Enumeroute(n).Ind(K(n),:)=Blink;
%     K(n)=K(n)+1;
% end
% K(n)=K(n)-1;
% KK(n)=K(n);
% 
% % Routes with n>1 links
% for n=2:25
%     K(n)=1;
%     for p=1:K(n-1)
%         for q=1:size(Conn(Enumeroute(n-1).Route(p,n)).Next,1)
%             if numel(intersect(Enumeroute(n-1).Route(p,:),Conn(Enumeroute(n-1).Route(p,n)).Next(q,1)))== 0
%                 Enumeroute(n).Route(K(n),:)=[Enumeroute(n-1).Route(p,:),Conn(Enumeroute(n-1).Route(p,n)).Next(q,1)];
%                 Enumeroute(n).TT(K(n),:)=[Enumeroute(n-1).TT(p,:),Conn(Enumeroute(n-1).Route(p,n)).TT(q,1)];
%                 Enumeroute(n).Cap(K(n),:)=[Enumeroute(n-1).Cap(p,:),Conn(Enumeroute(n-1).Route(p,n)).Cap(q,1)];
%                 Enumeroute(n).Ind(K(n),:)=[Enumeroute(n-1).Ind(p,:),Conn(Enumeroute(n-1).Route(p,n)).Ind(q,1)];
%                 K(n)=K(n)+1;
%             end
%         end
%     end
%     K(n)=K(n)-1;
%     KK(n)=sum(K(:));
% end

function[FinalRealDemand,thetat,RouteSet]=nlinknOD(LinkNode,Q,K,KK,Enumeroute)
nNode=size(K,2);
tau=20;
delta=0.25;
alpha=1+sqrt(log(2/delta)/2);
Theta=[2.5,-0.3];   % parameter of link function (binary logit model)
V=zeros(2,2);       % prepare variance matrix
X=zeros(1,2);
Y=zeros(1,1);
PreX=struct;
PreY=struct;     % prepare space for responses
Z=zeros(tau,nNode); % ???
PickedRoute=zeros(tau,1); % prepare space for initial route direction

% randomly choose route during learning period
for i=1:tau
    A(i,1)=randi([1,KK(nNode-1)]);  % choose random route ID among feasible route set
    for j=1:nNode-1
        if A(i,1)<=KK(j)
            %RouteLength=KK(j);  % ???
            PickedRoute(i,1)=j;  % identify length of chosen route
            break
        end
    end
    PickedRoute(i,2)=randi([1,K(PickedRoute(i,1))]);    % choose random route ID among feasible route set of chosen length
    
    % recall information about chosen route from enumerated route set
    for k=1:PickedRoute(i,1)+1
        PickNode(i,k)=Enumeroute(PickedRoute(i,1)).Route(PickedRoute(i,2),k);       % node ID
        PickQ(i,k)=Q(Enumeroute(PickedRoute(i,1)).Route(PickedRoute(i,2),k),1);     % potential demand for nodes in the sequence of route
        PickCumTT(i,k)=sum(Enumeroute(PickedRoute(i,1)).TT(PickedRoute(i,2),1:k));  % cumulative travel time
        PickCap(i,k)=Enumeroute(PickedRoute(i,1)).Cap(PickedRoute(i,2),k);          % capacity of included links
        PickInd(i,k)=Enumeroute(PickedRoute(i,1)).Ind(PickedRoute(i,2),k);          % index of included links
    end
end

% % simulate passengers' behaviors using binary logit model
% a=1;
% for i=1:tau
%     for j=1:PickedRoute(i,1)
%         for k=1:PickQ(i,j+1)
%             X(a,:)=[1,PickCumTT(i,j+1)];        % accumulate feature vectors
%             if rand<1/(1+exp(-Theta*X(a,:)'))   % compare random number with probability
%                 Y(a,1)=1;               % accepted 
%                 Z(i,PickNode(i,j+1))=Z(i,PickNode(i,j+1))+1;    % cumulative number of passengers to nodes
%             else
%                 Y(a,1)=0;               % rejected
%             end
%             V=V+X'*X;       % variance matrix
%             a=a+1;          % move on to next potential passenger
%         end
%     end
% end

% simulate passengers' behaviors using binary logit model (aggregated)
for i=1:tau
    for j=1:PickedRoute(i,1)
        a=1;
        for k=1:PickQ(i,j+1)
            PreX(PickNode(i,j+1)).X(i,:)=[1,PickCumTT(i,j+1)];        % accumulate feature vectors
            if rand<1/(1+exp(-Theta*PreX(PickNode(i,j+1)).X(i,:)'))   % compare random number with probability
                PreY(PickNode(i,j+1)).Y(i,a)=1;               % accepted 
                a=a+1;          % move on to next potential passenger
%                 Z(i,PickNode(i,j+1))=Z(i,PickNode(i,j+1))+1;    % cumulative number of passengers to nodes
            else
                PreY(PickNode(i,j+1)).Y(i,a)=0;               % rejected
            end
%             V=V+X'*X;       % variance matrix
        end
    end
end

% delete empty column
for i=1:nNode
    if size(PreY(i).Y,2)>0
    if sum(PreY(i).Y(:,size(PreY(i).Y,2)))==0
        PreY(i).Y(:,size(PreY(i).Y,2))=[];
    end
    end
end

% formulate X and Y for UCB-GLM, and V
b=0;
for i=1:nNode
    for j=1:size(PreX(i).X,1)
        if sum(PreX(i).X(j,1:2))>0
            for k=1:size(PreY(i).Y,2)
                b=b+1;
                X(b,:)=PreX(i).X(j,:);
                Y(b,1)=PreY(i).Y(j,k);
                V=V+X(b,:)'*X(b,:);
            end
        end
    end
end

% prepare space for archiving information about arms
ArmNode=zeros(KK(nNode-1),nNode);
ArmCumTT=zeros(KK(nNode-1),nNode);
ArmInd=zeros(KK(nNode-1),nNode);
ArmLink=zeros(KK(nNode-1),size(LinkNode,1));

% recall route information from struct and assign to matrix
c=1;
for i=1:nNode-1
    for j=1:K(i)
        routetime=0;
        for k=2:i+1
            if Enumeroute(i).TT(j,k)>0
                routetime=routetime+Enumeroute(i).TT(j,k);
                ArmNode(c,k)=Enumeroute(i).Route(j,k);
                ArmCumTT(c,Enumeroute(i).Route(j,k))=routetime;
                ArmInd(c,k)=Enumeroute(i).Ind(j,k);
                ArmLink(c,Enumeroute(i).Ind(j,k))=Enumeroute(i).Cap(j,k);
            end
        end
        c=c+1;
    end
end

budget=10;  % maximum number of routes

% prepare space for route set which consists of chosen routes
ArmUT=zeros(KK(nNode-1),nNode);
RouteSet=zeros(budget,1);
RouteSetInd=zeros(budget,nNode);
RouteSetNode=zeros(budget,nNode);
RouteSetCumTT=zeros(budget,nNode);
RouteSetUT=zeros(1,nNode);
RouteSetUTS=zeros(budget,nNode);
RouteSetODFlow=zeros(budget,nNode);
RouteSetLink=zeros(budget,size(ArmLink,2));
SumRouteSetLink=zeros(1,size(LinkNode,1));

RouteSetRealUT=zeros(1,nNode);
RouteSetRealUTS=zeros(budget,nNode);
RouteSetRealODFlow=zeros(budget,nNode);
RouteSetRealLink=zeros(budget,size(ArmLink,2));
SumRouteSetRealLink=zeros(1,size(LinkNode,1));

ResultLinkSet=zeros(KK(nNode-1),nNode);

% expand the route set by 1 route per trial
for t=tau+1:tau+budget
    % estimate parameters of logit model
    F=@(thetat)nlinknODfcn(thetat,X,Y,size(X,1));
    thetat=fsolve(F,zeros(1,2));
    theta(t,:)=thetat;
    
    
    tx=t-tau;
    TotalDemand(1:size(ArmCumTT,1),1)=0;
    TotalRealDemand(1:size(ArmCumTT,1),1)=0;
    maxd=-inf;
    maxRd=-inf;
    
    % evaluate the best route to add to set
    for m=1:size(ArmCumTT,1)
        if sum(RouteSet(:,1)==m)>0  % if the route already exists in the set
            continue
        else
            RouteSetODFlow(tx,:)=0;
            RouteSetNode(tx,:)=ArmNode(m,:);
            RouteSetLink(tx,:)=0;
            RouteSet(tx,1)=m;
            RouteSetUT(tx,:)=0;
            RouteSetCumTT(tx,:)=ArmCumTT(m,:);
            RouteSetInd(tx,:)=ArmInd(m,:);
            
            RouteSetRealODFlow(tx,:)=0;
            RouteSetRealUT(tx,:)=0;      
            RouteSetRealLink(tx,:)=0;
            
            % calculate the expected number of passengers by adding a route to the set
            d=0;
            for n=1:nNode
                for p=1:tx
                    if RouteSetCumTT(p,n)>0
                        RouteSetUT(p,n)=exp(thetat*[1;RouteSetCumTT(p,n)]); % calculate utility of route
                        RouteSetRealUT(p,n)=exp(Theta*[1;RouteSetCumTT(p,n)]); % calculate real utility of route
                    end
                end
                RouteSetUTS(tx,n)=sum(RouteSetUT(1:tx,n)); % sum of utilities for a node
                RouteSetRealUTS(tx,n)=sum(RouteSetRealUT(1:tx,n)); % sum of utilities for a node
                if RouteSetUTS(tx,n)~=0
                    for p=1:tx
                        RouteSetODFlow(p,n)=Q(n)*RouteSetUT(p,n)/(1+RouteSetUTS(tx,n)); % calculate portion of potential demand 
                        RouteSetRealODFlow(p,n)=Q(n)*RouteSetRealUT(p,n)/(1+RouteSetRealUTS(tx,n)); % calculate real portion of potential demand 
                    end
                end
            end

            for p=1:tx-1
                n=2;
%                 RouteSetLink(p,RouteSetInd(p,n))=sum(RouteSetODFlow(p,:));
                RouteSetLink(p,RouteSetInd(p,n))=sum(RouteSetRealODFlow(p,:));
                RouteSetRealLink(p,RouteSetInd(p,n))=sum(RouteSetRealODFlow(p,:));
                n=3;
                while RouteSetInd(p,n)>0
                    RouteSetLink(p,RouteSetInd(p,n))=RouteSetLink(p,RouteSetInd(p,n-1))-RouteSetRealODFlow(p,RouteSetNode(p,n-1));
%                     RouteSetLink(p,RouteSetInd(p,n))=RouteSetLink(p,RouteSetInd(p,n-1))-RouteSetODFlow(p,RouteSetNode(p,n-1));
                    RouteSetRealLink(p,RouteSetInd(p,n))=RouteSetRealLink(p,RouteSetInd(p,n-1))-RouteSetRealODFlow(p,RouteSetNode(p,n-1));
                    n=n+1;
                    if n==nNode
                        break
                    end
                end
            end
            
            for p=tx
                n=2;
                RouteSetLink(p,RouteSetInd(p,n))=sum(RouteSetODFlow(p,:));
                RouteSetRealLink(p,RouteSetInd(p,n))=sum(RouteSetRealODFlow(p,:));
                n=3;
                while RouteSetInd(p,n)>0
                    RouteSetLink(p,RouteSetInd(p,n))=RouteSetLink(p,RouteSetInd(p,n-1))-RouteSetODFlow(p,RouteSetNode(p,n-1));
                    RouteSetRealLink(p,RouteSetInd(p,n))=RouteSetRealLink(p,RouteSetInd(p,n-1))-RouteSetRealODFlow(p,RouteSetNode(p,n-1));
                    n=n+1;
                    if n==nNode
                        break
                    end
                end
            end

            for r=1:size(LinkNode,1)
                SumRouteSetLink(1,r)=sum(RouteSetLink(:,r));
                SumRouteSetRealLink(1,r)=sum(RouteSetRealLink(:,r));
                if SumRouteSetLink(1,r)<LinkNode(r,4)
                    TotalDemand(m,1)=sum(sum(RouteSetODFlow(1:tx,:)));
                    TotalRealDemand(m,1)=sum(sum(RouteSetRealODFlow(1:tx,:)));
                else
                    TotalDemand(m,1)=-inf;
                    TotalRealDemand(m,1)=-inf;
                    break
                end
            end
            if TotalDemand(m,1)>maxd
                maxd=TotalDemand(m,1);
                BB=RouteSetNode(tx,:);
                CC=RouteSetLink(tx,:);
                DD=RouteSetODFlow(tx,:);
                EE=RouteSetUTS(tx,:);
                FF=RouteSetCumTT(tx,:);
                GG=RouteSetUT(tx,:);
                
                LL=RouteSetRealUT(tx,:);
                HH=RouteSetRealUTS(tx,:);
                JJ=RouteSetRealODFlow(tx,:);
                MM=RouteSetRealLink(tx,:);
                maxRd=TotalRealDemand(m,1);
            end
        end
        if maxd==-inf
            break
        end
    end
    if maxd==-inf
        break
    end
    [~,I]=max(TotalDemand);
    RouteSet(tx,1)=I;
    RouteSetNode(tx,:)=BB;
%     RouteSetLink(tx,:)=CC;
    RouteSetLink(tx,:)=MM;
    RouteSetUT(tx,:)=GG;
    RouteSetUTS(tx,:)=EE;
    RouteSetODFlow(tx,:)=DD;
    RouteSetCumTT(tx,:)=FF;
    RouteSetInd(tx,:)=ArmInd(I,:);
    
    RouteSetRealODFlow(tx,:)=JJ;
    RouteSetRealUTS(tx,:)=HH;
    RouteSetRealUT(tx,:)=LL;
    RouteSetRealLink(tx,:)=MM;

    % simulate passengers' behaviors using binary logit model (aggregated)
    
%         for j=1:PickedRoute(i,1)
%     for i=1:tx
        for j=1:sum(RouteSetNode(tx,:)>0)
            a=b+1;
            for k=1:Q(RouteSetNode(tx,j+1),1)
%                 b=b+1;
                ArmValue=[];
                PreX(RouteSetNode(tx,j+1)).X(t,:)=[1,RouteSetCumTT(tx,RouteSetNode(tx,j+1))];        % accumulate feature vectors
                if size(PreX(RouteSetNode(tx,j+1)).X(:,1),1)>t-1
                    for h=1:size(PreX(RouteSetNode(tx,j+1)).X(:,1),1)-tau
                        if PreX(RouteSetNode(tx,j+1)).X(tau+h,1)==1
                            ArmValue(h,1)=thetat*PreX(RouteSetNode(tx,j+1)).X(t,:)'+alpha*sqrt(PreX(RouteSetNode(tx,j+1)).X(t,:)*V^-1*PreX(RouteSetNode(tx,j+1)).X(t,:)');
                            ArmValue(h,2)=RouteSetNode(tx,j+1);
                        else
                            ArmValue(h,1)=-inf;
                            ArmValue(h,2)=0;
                        end
                    end
                    [maxArmValue,maxArm]=max(ArmValue(:,1));
                end
                
                
                if rand<1/(1+1/RouteSetRealUT(tx,ArmValue(maxArm,2)))   % compare random number with probability
                    PreY(RouteSetNode(tx,j+1)).Y(t,a)=1;               % accepted 
                    Y(a,1)=1;
                    X(a,:)=PreX(ArmValue(maxArm,2)).X(tau+maxArm,:);
                    V=V+X(a,:)'*X(a,:);
                    a=a+1;          % move on to next potential passenger
                else
                    PreY(RouteSetNode(tx,j+1)).Y(t,a)=0;               % rejected
                end
            end
        end 
%         for j=1:sum(RouteSetNode(tx,:)>0)
%             a=1;
%             for k=1:Q(RouteSetNode(tx,j+1),1)
%                 PreX(RouteSetNode(tx,j+1)).X(t,:)=[1,RouteSetCumTT(tx,RouteSetNode(tx,j+1))];        % accumulate feature vectors
%                 if rand<1/(1+1/RouteSetRealUTS(tx,RouteSetNode(tx,j+1)))   % compare random number with probability
%                     PreY(RouteSetNode(tx,j+1)).Y(t,a)=1;               % accepted 
%                     a=a+1;          % move on to next potential passenger
%     %                 Z(i,PickNode(i,j+1))=Z(i,PickNode(i,j+1))+1;    % cumulative number of passengers to nodes
%                 else
%                     PreY(RouteSetNode(tx,j+1)).Y(t,a)=0;               % rejected
%                 end
%             end
%         end
%     end
%     
% 
%     % delete empty column
%     for i=1:nNode
%         if size(PreY(i).Y,2)>0
%             if sum(PreY(i).Y(:,size(PreY(i).Y,2)))==0
%                 PreY(i).Y(:,size(PreY(i).Y,2))=[];
%             end
%         end
%     end
%     % formulate 
%     for i=1:nNode
%         SumRouteSetRealUTS(tx,i)=sum(RouteSetRealUTS(:,i));
%     end
% %     formulate X and Y for UCB-GLM
%     for i=1:nNode
%         if size(PreX(i).X,1)==t && sum(PreX(i).X(t,1:2))>0
%             for k=1:size(PreY(i).Y,2)
%                 b=b+1;
%                 p=rand;
%                 q=1;
%                 while p>SumRouteSetRealUTS(q,i)/SumRouteSetRealUTS(tx,i) && q<=tx
%                     q=q+1;
%                 end
%                 X(b,:)=PreX(i).X(tau+q,:);                    
%                 Y(b,1)=PreY(i).Y(t,k);
%                 V=V+X(b,:)'*X(b,:);
%             end
%         end
%     end   
end

FinalDemand=sum(sum(RouteSetODFlow(:,:)))
% FinalRealDemand=sum(sum(RouteSetRealODFlow(:,:)))

XXX=zeros(budget,nNode);
for i=1:nNode
    for j=1:budget
        XXX(j,i)=Q(i)*RouteSetRealUT(j,i)/(1+RouteSetRealUTS(budget,i));
    end
end
FinalRealDemand=sum(sum(XXX(:,:)))

end