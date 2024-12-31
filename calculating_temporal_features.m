% load preprocessed fMRI time series 

nscan = size(FC,1); 
nroi  = size(FC,2);
nsubj = size(FC,3);

%ide are functional connectivities between amygdala, nucleus accumbens, and other
%regions

nedge = length(ide); 
% calculate functional connectivity (FC)
RS = zeros(nsubj,nedge);
for s=1:nsubj
    R1 = corr(FC(:,:,s));       % individual functional network
    RS(s,:) = R1(ide);          % lower triangular elements
end
ZS = atanh(RS) * sqrt(nscan-3); % ZS ~ N(0,1) : Fisher's r-to-z

% normalize each edge to be N(0,1) 
M  = ZS;
EM = repmat(mean(M),size(M,1),1);
SM = repmat(std(M), size(M,1),1);
Z  = (M - EM) ./ SM;

[ncomp, logp] = laplace_pca(Z); 

% PCA and ICA on networks

ncomp = 43;%%%numbers using laplace (first analysis)

[X,dwM,L,V,whiteM,V_orig,L_orig] = rsn_pca(Z',ncomp);  
[iq,A_ica,W_ica,S_ica,sR_ica]    = icasso(double(X)',30,'g','tanh', ...
                                  'lastEig',ncomp,'vis','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% z-normalization for thresholding
GS=S_ica'; % nedge x ncomp
MS=repmat(mean(GS),size(GS,1),1);
VS=repmat(std(GS),size(GS,1),1);
GS= (GS - MS) ./ VS; % ~ N(0,1)

thr = 3;             % edge-level threshold height (z-score scale)
idcomp=find(iq>0.80); % iq = stability 
S = GS(:,idcomp); 
S(abs(S)<=thr)=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for sequence making and visualization
wsize = 42; % 0.72 x 42 = 30.24 secs
Y = zeros(nscan-wsize,size(GS,2),nsubj); % time x comp x subj
for s=1:nsubj
    for t=1:nscan-wsize
        Yw  = FC(t:t+wsize,:,s);
        Rw  = corr(Yw);
        Zw1 = atanh(Rw) * sqrt(wsize-3);
        ZW  = Zw1(ide);  % FC for the current window
        
        idinf = find(isinf(ZW)); %% if correlation 1? inf -> should be changed for max value
        ZW(idinf) = 0; 
        ZW(idinf) = max(ZW); 
        
        ZW  = (ZW - mean(ZW)) / std(ZW); % intensity normalization

        % linear projection on IC basis space
        Y(t,:,s)=pinv(GS)*ZW; % OLS approach
        if rem(t,100)==0
            disp([num2str(t) '-th window was done ...']);
        end
    end
    disp([num2str(s) '-th subject was done ...']);    
end

%%selected IC maps
%Y_pos = Y(:,[choose maps],:); 
%Y_neg = Y(:,[choose maps],:);

%change the negative signal
Y_neg = -Y_neg;


for s=1:1080
    Y_c(:,:,s) = horzcat(Y_pos(:,:,s),Y_neg(:,:,s));
end

%%%%%%%%%%%%%%%%%%%%for calculating variability
%%for variance

index_var = zeros(1080,9);
for s=1:1080
    variance_values  = var(Y_c(:,:,s));
    index_var(s,:) = variance_values; 
    disp([num2str(s) '-th subject was done ...']);    
end


%%for rms

index_rms = zeros(1080,9);
for s=1:1080
    rms_values  = sqrt(mean(Y_c(:,:,s).^2)); 
    index_rms(s,:) = rms_values; 
    disp([num2str(s) '-th subject was done ...']);    
end







