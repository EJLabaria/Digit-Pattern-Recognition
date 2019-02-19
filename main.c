% Singular Value Decomposition Algorithm
 
%  (1) Sort into a 3 x 3 matrix/array
%  (2) Once each digit is sorted turn each digit into a 2 x 2
%  (3) Run SVD on each digit syntax: [Ux V D] = svd(tensorx)
%         where x is a specific digit for specific tensor 
%   Run residuals algorithm
%   Classify unkown digit
 
 
A = azip;
D = dzip;
test = testzip;
 
%%% INITIALIZATION   %%%
sorting = zeros(256, 10 , 1707);  
ecount = zeros(10,1);
%tensor = zeros(256,1707);
 
% Sorting loop
for i = 1:1707
    for j = 0:9;
        
        if D(i) == j;   
            ecount(j+1)= ecount(j+1) +1; 
           sorting(:, j+1 , ecount(j+1)) = A(:,i);
           
        end % end
    end % end
end % end 
 
 
 
% NOW CONVERT EACH DIGIT TO A 2X2 INSTEAD OF A 3X3 BASICALLY SPERATE EACH
% DIGIT INTO ITS OWN MATRIX INSTEAD OF A BIG BOX OF SLICES
 
 
% k = 18 = 94.3199
 
k = 17;
 
for z= 1:10
    if z == 1
        for i = 1:ecount(z)
            prac0(:,i) = sorting(:,z,i);
        end
        
        [U0 S V] = svds(prac0,k);
       
    
        elseif z == 2
        for i = 1:ecount(z)
            prac1(:,i) = sorting(:,z,i);
        end
        
        [U1 S V] = svds(prac1,k);
        
        elseif z == 3
        for i = 1:ecount(z)
            prac2(:,i) = sorting(:,z,i);
        end
        
        [U2 S V] = svds(prac2,k);
        
       elseif z == 4
        for i = 1:ecount(z)
            prac3(:,i) = sorting(:,z,i);
        end 
        [U3 S V] = svds(prac3,k);
      
       elseif z == 5
        for i = 1:ecount(z)
            prac4(:,i) = sorting(:,z,i);
        end 
        [U4 S V] = svds(prac4,k);
 
        elseif z == 6
        for i = 1:ecount(z)
            prac5(:,i) = sorting(:,z,i);
        end 
        
        [U5 S V] = svds(prac5,k);
        
        elseif z == 7
  
            for i = 1:ecount(z)
            prac6(:,i) = sorting(:,z,i);
            end 
            [U6 S V] = svds(prac6,k);
        
    elseif z == 8
        for i = 1:ecount(z)
            prac7(:,i) = sorting(:,z,i);
        end 
    [U7 S V] = svds(prac7,k);
    
        elseif z == 9
        for i = 1:ecount(z)
            prac8(:,i) = sorting(:,z,i);
        end 
        [U8 S V] = svds(prac8,k);
   
    elseif z ==10
            for i = 1:ecount(z)
             prac9(:,i) = sorting(:,z,i);
            end
        [U9 S V] = svds(prac9,k);    
            
    end % end of IF loop
end % END of OUTER Loop    
    
 
%%% Compute U*U' to minimuze calculation/time
 
U0t = U0*U0';
U1t = U1*U1';
U2t = U2*U2';
U3t = U3*U3';
U4t = U4*U4';
U4t = U4*U4';
U5t = U5*U5';
U6t = U6*U6';
U7t = U7*U7';
U8t = U8*U8';
U9t = U9*U9';
 
 
% Compute residuals between test vectors 
res = zeros(10,size(test,2));
results2 = zeros(1,size(test,2));
suc2 = 0;
for i = 1:size(testzip,2)
    % must compare each test image to all ten basis vectors
    for j = 1:10
        if j == 1
            res(j,i) = norm((eye(size(U0t,1)) - U0t) * test(:,i), 2);
        elseif j == 2
            res(j,i) = norm((eye(size(U1t,1)) - U1t) * test(:,i), 2);
        elseif j == 3
            res(j,i) = norm((eye(size(U2t,1)) - U2t) * test(:,i), 2);
        elseif j == 4
            res(j,i) = norm((eye(size(U3t,1)) - U3t) * test(:,i), 2);
        elseif j == 5
            res(j,i) = norm((eye(size(U4t,1)) - U4t) * test(:,i), 2);
        elseif j == 6
            res(j,i) = norm((eye(size(U5t,1)) - U5t) * test(:,i), 2);
        elseif j == 7
            res(j,i) = norm((eye(size(U6t,1)) - U6t) * test(:,i), 2);
        elseif j == 8
            res(j,i) = norm((eye(size(U7t,1)) - U7t) * test(:,i), 2);
        elseif j == 9
            res(j,i) = norm((eye(size(U8t,1)) - U8t) * test(:,i), 2);
        elseif j == 10
            res(j,i) = norm((eye(size(U9t,1)) - U9t) * test(:,i), 2);
        end
    end
    
    min = res(1,i);
    results2(i) = 0;
    for k = 2:10
       if res(k,i) < min
           min = res(k,i);
           results2(i) = (k - 1);
       end
    end
    
    if results2(i) == dtest(i)
       suc2 = suc2 + 1; 
    end
end
 
successRate2 = 100 * suc2 / size(testzip,2);
disp('Rate of success of');
disp('SVD classification:');
disp(successRate2);
