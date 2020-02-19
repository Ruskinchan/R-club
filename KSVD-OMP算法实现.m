function [A]=OMPerr(D,X,errorGoal)
[n,P]=size(X);%n=64  P= 62001=249*249
[n,K]=size(D);%n=64 K=256
E2 = errorGoal^2*n;
maxNumCoef = n/2;%%%%%%32
A = sparse(size(D,2),size(X,2));
for k=1:1:P,
    a=[];
    x=X(:,k);
    residual=x;
    indx = [];
    a = [];
    currResNorm2 = sum(residual.^2);
    j = 0;
    while currResNorm2>E2 & j < maxNumCoef,
        j = j+1;
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos; 
        a=pinv(D(:,indx(1:j)))*x;%j*64  *64*1=j*1    
        residual=x-D(:,indx(1:j))*a;
        currResNorm2 = sum(residual.^2);
   end;
   if (length(indx)>0)
       A(indx,k)=a;
   end
end;
return;
