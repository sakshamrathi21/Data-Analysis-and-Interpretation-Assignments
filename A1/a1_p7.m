clear; clc;
rng(42);

a = [1 2 4 6 7 2];
b = [1 2 4 6 7 2 4];
disp('mean test:');
disp(mean(b));
disp(UpdateMean(mean(a), 4, length(a)));

disp('std test:');
disp(std(b));
disp(UpdateStd(mean(a), std(a), UpdateMean(mean(a), 4, length(a)), 4, length(a)));

n = 10;
m = 11;
c = rand(n, 1);
d = zeros(m, 1);

for i = 1:n
    d(i) = c(i);
end

new_elem = rand();
d(m) = new_elem;

disp('median test:');
disp(median(d));

sorted_c = sort(c);
disp(UpdateMedian(new_elem, sorted_c, length(sorted_c)));

function newMean = UpdateMean(OldMean, NewDataValue, n)
    newMean = (n * OldMean + NewDataValue) / (n + 1);
end

function newMedian = UpdateMedian(NewDataValue, A, n)
    if(mod(n, 2) == 0)
        % n is even
        mid_left = A(n/2);
        mid_right = A(n/2 + 1);
        
        if(NewDataValue < mid_left)
            newMedian = mid_left;
        elseif(NewDataValue > mid_right)
            newMedian = mid_right;
        else
            newMedian = NewDataValue;
        end
    else
        % n is odd
        mid = A((n+1)/2);
        mid_left = A((n-1)/2);
        mid_right = A((n+1)/2 + 1);

        if(NewDataValue < mid_left)
            newMedian = median(mid_left, mid);
        elseif(NewDataValue > mid_right)
            newMedian = median(mid, mid_right);
        else
            newMedian = median(NewDataValue, mid);
        end
    end
end

function newStd = UpdateStd (OldMean, OldStd, NewMean, NewDataValue, n)
    variance = ( (n-1) * OldStd^2 + NewDataValue^2 + n * OldMean^2 - (n+1) * NewMean^2 ) / n;
    newStd = sqrt(variance);
end
