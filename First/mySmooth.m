function res = mySmooth(x, m)
    n = size(x, 2);
    tx = zeros(m, n);
    tx1 = [tx; x];
    tx2 = [x; tx];
    res = (tx1 + tx2) / m;
    
end