function H = entropy(P)


%P = sum(Pxr, dim); % sum over dim
P = P(find(P)); % removes zeros
P = P(find(isnan(P)==0)); % removes NaNs
P = P(find(isinf(P)==0)); % removes Infs
H = -sum(P.*log2(P)); % in bits