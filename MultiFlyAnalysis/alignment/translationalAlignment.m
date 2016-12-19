function [X,Y] = translationalAlignment(image1,image2)

    s = size(image1);
    Q = real(ifft2(fft2(image1).*conj(fft2(image2))));
    [ii,jj] = find(Q == max(Q(:)));
    
    xIdx = mod0((ii-1):(ii+1),s(1));
    yIdx = mod0((jj-1):(jj+1),s(2));
    
    Z = Q(xIdx,yIdx);
    
    A = [0 0 0 0 0 1;
         1 0 0 1 0 1;
         0 1 0 0 1 1;
         1 0 0 -1 0 1;
         0 1 0 0 -1 1;
         1 1 1 1  1 1];
     
     Aends = [1 1 1 1 1 1;
         1 1 -1 1 -1 1;
         1 1 1 -1 -1 1;
         1 1 -1 -1 1 1];
     
     Zends = [Z(3,3);-Z(3,1);Z(1,1);-Z(1,3)];
     
     vals = zeros(6,4);
     for i=1:4
         A(6,:) = Aends(i,:);
         vals(:,i) = A \ [Z(2,2);Z(3,2);Z(2,3);Z(1,2);Z(2,1);Zends(i)];
     end
     
     x = mean(vals,2);
     
     a = x(1);
     b = x(2);
     c = x(3);
     d = x(4);
     e = x(5);
     
     X = jj - (2*b + e)/c;
     Y = ii - (2*a + d)/c;
     

     if X > s(1)/2
         X = X - s(1) - 1;
     else
         X = X - 1;
     end
     
     if Y > s(2)/2
         Y = Y - s(2) - 1;
     else
         Y = Y - 1;
     end

     
     if length(X) > 1
         X = X(1);
     else
         if isempty(X)
             X = 0;
         end
     end
     
     
     if length(Y) > 1
         Y = Y(1);
     else
         if isempty(Y)
             Y = 0;
         end
     end
     