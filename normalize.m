function [A] = normalize(A)

[~,n]=size(A);

for i=1:n
    a=sqrt(A(:,i)'*A(:,i));
    A(:,i)=A(:,i)./a;
end

end

