function val = sCore(n,k)
%     d = ((n-k)/k + k-1);
%     val = (n-k)*d + (k*(k-1)/2)*d^2;
    
    d1 = ((n-k)/k + k-1);
    d2 = ((n-k)/k + 1);
    val = ((n-k)/k)*d1 + (k-1)*((n-k)/k)*d2 + (k-1)*d1*d2;