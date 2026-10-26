function H = gyro_sampen(x, maxPoints, rScale)
% Exact SampEn on a uniformly sub-sampled vector (m=2, r=rScale*std).
% Uses pair counts excluding self-matches; maxPoints bounds O(N^2).
x=double(x(:));
if numel(x)>maxPoints
    ix=round(linspace(1,numel(x),maxPoints)); x=x(ix);
end
n=numel(x); r=rScale*std(x);
if n<12 || r<=eps
    H=0; return
end
B=0; A=0;
% Only vector pairs with enough points for both m and m+1.
for i=1:n-3
    j=(i+1):(n-2);
    close2=(abs(x(i)-x(j))<=r) & (abs(x(i+1)-x(j+1))<=r);
    B=B+sum(close2);
    A=A+sum(close2 & (abs(x(i+2)-x(j+2))<=r));
end
if B==0
    H=Inf;
elseif A==0
    H=Inf;
else
    H=-log(A/B);
end
end
