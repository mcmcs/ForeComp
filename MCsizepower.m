%% Lazarus, Lewis, and Stock (2020), "The Size-Power Tradeoff in HAR Inference"
% This function runs Monte Carlo simulations for a given DGP and calculates
% size and the maximum power loss using each of the 5 (six including KVB)
% estimators considered in the paper.
% Inputs:
% alpha: ARMA parameter vector for DGP
% T: sample size
% DGP: @name_of_dgp, calls DGP
% draws: number of Monte Carlo samples
% regs: dimension of parameter vector of interest. For location model, this
% equals m. For stochastic regressor model, it is the number of stochastic
% regressors.
% Outputs:
% MCsize: vector of Monte Carlo sizes for each test
% maxPloss: vector of Monte Carlo worst-case power losses for each test

function [MCsize,maxPloss]=MCsizepower(alpha,T,DGP,draws,regs,m)

% Compute truncation parameters etc.
dof=[8,16,32,64];
BQS=ceil(dof/(5/3));
SNW=ceil(3*T./(2*dof));
Bfourier=dof/2;

% Construct kernels and basis functions
kkvb=kvbkernel(T); %Time domain
Wqs=[QSfreq(T,BQS(1));QSfreq(T,BQS(2));QSfreq(T,BQS(3));QSfreq(T,BQS(4))]; %Frequency domain
knw=[nwkernel(T,SNW(1));nwkernel(T,SNW(2));nwkernel(T,SNW(3));nwkernel(T,SNW(4))]; %Time domain
fbasis=wpbasis(T,max(Bfourier));
cbasis=cosbasis(T,max(dof));
ssbasis=[imbasis(T,dof(1));imbasis(T,dof(2));imbasis(T,dof(3));imbasis(T,dof(4))];

% Construct Storage
betadev=NaN(draws,regs);
omegastoreQS=NaN(draws,4,regs,regs);
omegastoreNW=NaN(draws,4,regs,regs);
omegastoreEWP=NaN(draws,4,regs,regs);
omegastoreCOS=NaN(draws,4,regs,regs);
omegastoreSS=NaN(draws,4,regs,regs);
omegastoreKVB=NaN(draws,1,regs,regs);
ACFind=reshape(reshape(1:m^2,m,m)',m^2,1);
% Seed
rng(715);

for i=1:draws
    
    % Generate data
    [zeta,b,~,data,~]=DGP(T,alpha,m);
    % Store betadev
    betadev(i,:)=b;
    
    % Compute ACF
    ACF=NaN(T-1,regs^2);
    for j=1:T
        ACF(j,:)=reshape(((zeta(j:end,:))'*(zeta(1:end-j+1,:))),regs^2,1)/T; %Compute ACF at lag j without substracting mean - as per Jim's notes
    end
    if regs==1
        ACF=[flipud(ACF(2:end,:));ACF];
        
    else
        ACF=[flipud(ACF(2:end,:));ACF(:,ACFind)]; % Duplicate for -j
    end
    XX=data'*data/T; %Compute Sig_XX
    
    % Calculate LRVs
    omegastoreKVB(i,1,:,:)=XX^-1*reshape([fliplr(kkvb),1,kkvb]*ACF,regs,regs)*XX^-1;
    for d=1:length(dof)
        if m==1
            FT=fft(zeta);
            omega=Wqs(d,1:BQS(d))*(FT(2:BQS(d)+1,:).*conj(FT(2:BQS(d)+1,:))/T);
        else
            % Accounting for the non-unitary Izz ordinates
            omega=zeros(m); % will collect the summation of weights*Izz
            FT=fft(zeta); % compute FT of zeta, dz
            for l=1:BQS(d) % loop through the ordinates of the periodogram, up to J=BQS(d)
                Izz=FT(l+1,:)'*ctranspose(FT(l+1,:)')/T; % Compute the l+1th ordinate, Izz=dz(w)*conj(dz(w))'/T
                omega=omega+Wqs(d,l)*real(Izz); % add the lth weight times this ordinate to the sum
            end
        end
        omegastoreQS(i,d,:,:)=(XX^-1*omega*XX^-1); % Compute the variance using X'X and omega
        omegastoreNW(i,d,:,:)=XX^-1*reshape([fliplr(knw(d,1:(length(ACF)-1)/2)),1,knw(d,1:(length(ACF)-1)/2)]*ACF,regs,regs)*XX^-1;
        h=fbasis(1:dof(d),1:size(zeta,1))*zeta/T^.5;
        omegastoreEWP(i,d,:,:)=(XX^-1*1/dof(d)*(h'*h)*XX^-1);
        h=cbasis(1:dof(d),1:size(zeta,1))*zeta/T^.5;
        omegastoreCOS(i,d,:,:)=(XX^-1*1/dof(d)*(h'*h)*XX^-1);
        h=ssbasis(sum(dof(1:d-1))+1:sum(dof(1:d)),1:size(zeta,1))*zeta/T^.5;
        omegastoreSS(i,d,:,:)=(XX^-1*1/dof(d)*(h'*h)*XX^-1);
    end
    
end
%% Size tests
if m==1
    % Construct t stats
    tQS=T^.5*repmat(betadev,1,length(dof))./squeeze(omegastoreQS(:,:,1,1).^.5);
    tNW=T^.5*repmat(betadev,1,length(dof))./squeeze(omegastoreNW(:,:,1,1).^.5);
    tEWP=T^.5*repmat(betadev,1,length(dof))./squeeze(omegastoreEWP(:,:,1,1).^.5);
    tCOS=T^.5*repmat(betadev,1,length(dof))./squeeze(omegastoreCOS(:,:,1,1).^.5);
    tSS=T^.5*repmat(betadev,1,length(dof))./squeeze(omegastoreSS(:,:,1,1).^.5);
    tKVB=T^.5*betadev./squeeze(omegastoreKVB(:,:,1,1).^.5);
    
    
    % Critical values
    QScrits=[-2.3386,2.3375;-2.1203,2.1210;-2.0345,2.0338;-1.9963,1.9977,]'; %% Externally calculated, see readme
    QSmat1=repmat(QScrits(1,:),draws,1);
    QSmat2=repmat(QScrits(2,:),draws,1);
    NWcrits=[-2.4918,2.4924;-2.2200,2.2135;-2.0896,2.0797;-2.0219,2.0233]'; %% Externally calculated, see readme
    NWmat1=repmat(NWcrits(1,:),draws,1);
    NWmat2=repmat(NWcrits(2,:),draws,1);
    KVBcrits=[-4.771,4.771]'; % see KVB note
    KVBmat1=repmat(KVBcrits(1,:),draws,1);
    KVBmat2=repmat(KVBcrits(2,:),draws,1);
    basisfncrits=[tinv(.025,dof'),tinv(.975,dof')]'; % exact for fixed-b
    basismat1=repmat(basisfncrits(1,:),draws,1);
    basismat2=repmat(basisfncrits(2,:),draws,1);
    
    % Rejection Rates
    QSrej=mean((1-(QSmat1<tQS & tQS<QSmat2)));
    NWrej=mean((1-(NWmat1<tNW & tNW<NWmat2)));
    EWPrej=mean((1-(basismat1<tEWP & tEWP<basismat2)));
    COSrej=mean((1-(basismat1<tCOS & tCOS<basismat2)));
    SSrej=mean((1-(basismat1<tSS & tSS<basismat2)));
    KVBrej=mean((1-(KVBmat1<tKVB & tKVB<KVBmat2)));
    MCsize=[QSrej,NWrej,EWPrej,COSrej,SSrej,KVBrej];
    
    
    % Quantiles (and store in reference matrices for use in SAP
    % calculations)
    QSquantiles=quantile(tQS,[.025,.975]);
    QSmat1=repmat(QSquantiles(1,:),draws,1);
    QSmat2=repmat(QSquantiles(2,:),draws,1);
    NWquantiles=quantile(tNW,[.025,.975]);
    NWmat1=repmat(NWquantiles(1,:),draws,1);
    NWmat2=repmat(NWquantiles(2,:),draws,1);
    EWPquantiles=quantile(tEWP,[.025,.975]);
    EWPmat1=repmat(EWPquantiles(1,:),draws,1);
    EWPmat2=repmat(EWPquantiles(2,:),draws,1);
    COSquantiles=quantile(tCOS,[.025,.975]);
    COSmat1=repmat(COSquantiles(1,:),draws,1);
    COSmat2=repmat(COSquantiles(2,:),draws,1);
    SSquantiles=quantile(tSS,[.025,.975]);
    SSmat1=repmat(SSquantiles(1,:),draws,1);
    SSmat2=repmat(SSquantiles(2,:),draws,1);
    KVBquantiles=quantile(tKVB,[.025,.975]);
    KVBmat1=repmat(KVBquantiles(1),draws,1);
    KVBmat2=repmat(KVBquantiles(2),draws,1);
    %% Construct power curves
    
    % Generate alternatives
    delrange=0:.25:4;
    if isequal(DGP,@ARmean)
        betas=delrange/(1-alpha)/sqrt(T);
    elseif isequal(DGP,@stochasAR)
        rho=alpha^2;
        betas=delrange*sqrt((1+rho)/(1-rho))/sqrt(T);
    elseif isequal(DGP,@ARMA21)
        betas=delrange*sqrt(1.7668)/sqrt(T); % Spectral density from Figure 6 etc.
    end
    
    % Storage
    QSpow=NaN(length(delrange),4);
    NWpow=NaN(length(delrange),4);
    EWPpow=NaN(length(delrange),4);
    COSpow=NaN(length(delrange),4);
    SSpow=NaN(length(delrange),4);
    KVBpow=NaN(length(delrange),1);
    
    for j=1:length(delrange);
        betatest=betadev+repmat(betas(j),draws,1); % Augment estimated beta deviation from run with beta_alt to calculate power under alterntive
        % Test stats
        tQS=T^.5*repmat(betatest,1,length(dof))./squeeze(omegastoreQS(:,:,1,1).^.5);
        tNW=T^.5*repmat(betatest,1,length(dof))./squeeze(omegastoreNW(:,:,1,1).^.5);
        tEWP=T^.5*repmat(betatest,1,length(dof))./squeeze(omegastoreEWP(:,:,1,1).^.5);
        tCOS=T^.5*repmat(betatest,1,length(dof))./squeeze(omegastoreCOS(:,:,1,1).^.5);
        tSS=T^.5*repmat(betatest,1,length(dof))./squeeze(omegastoreSS(:,:,1,1).^.5);
        tKVB=T^.5*betatest./squeeze(omegastoreKVB(:,:,:,1,1).^.5);
        % Rejection rates
        QSpow(j,:)=mean((1-(QSmat1<tQS & tQS<QSmat2)));
        NWpow(j,:)=mean((1-(NWmat1<tNW & tNW<NWmat2)));
        EWPpow(j,:)=mean((1-(EWPmat1<tEWP & tEWP<EWPmat2)));
        COSpow(j,:)=mean((1-(COSmat1<tCOS & tCOS<COSmat2)));
        SSpow(j,:)=mean((1-(SSmat1<tSS & tSS<SSmat2)));
        KVBpow(j,:)=mean((1-(KVBmat1<tKVB & tKVB<KVBmat2)));
    end
    LGpow=NaN(size(SSpow));
    powsample=[QSpow,NWpow,EWPpow,COSpow,SSpow,KVBpow];
    
    %% Find worst case power
    grid=2:.05:4;
    samp=delrange;
    pow_gau = normcdf(-1.96+grid)+normcdf(-1.96-grid); % Oracle power
    powinterp=interp1(samp, powsample,grid,'linear'); % Interpolate power curves based on simulations
    [maxPloss]=max(repmat(pow_gau',1,length(MCsize))-powinterp); % Find maximum power loss
    
else % for m=2
    %% Size tests
    
    factor=(dof-m+1)./dof; % See Sun (2013)
    
    % Storage
    FQS=zeros(draws,length(dof));
    FNW=zeros(draws,length(dof));
    FEWP=zeros(draws,length(dof));
    FCOS=zeros(draws,length(dof));
    FSS=zeros(draws,length(dof));
    FKVB=zeros(draws,1);
    
    % Test statistics
    for i=1:draws
        for d=1:length(dof)
            FQS(i,d)=T*betadev(i,:)*squeeze(omegastoreQS(i,d,:,:))^-1*betadev(i,:)'/m;
            FNW(i,d)=T*betadev(i,:)*squeeze(omegastoreNW(i,d,:,:))^-1*betadev(i,:)'/m;
            FEWP(i,d)=factor(d).*T*betadev(i,:)*squeeze(omegastoreEWP(i,d,:,:))^-1*betadev(i,:)'/m;
            FCOS(i,d)=factor(d).*T*betadev(i,:)*squeeze(omegastoreCOS(i,d,:,:))^-1*betadev(i,:)'/m;
            FSS(i,d)=factor(d).*T*betadev(i,:)*squeeze(omegastoreSS(i,d,:,:))^-1*betadev(i,:)'/m;
        end
        FKVB(i)=T*betadev(i,:)*squeeze(omegastoreKVB(i,1,:,:))^-1*betadev(i,:)'/m;
    end
    
    % Critical values
    if m==2
        QScrits=[5.6501    3.9264    3.3978    3.1943]; % Externally calculated, see readme
        QSmat1=repmat(QScrits(1,:),draws,1);
        NWcrits=[5.9790,4.2755,3.5859,3.2770]; % Externally calculated, see readme
        NWmat1=repmat(NWcrits(1,:),draws,1);
        KVBcrits=51.41/2; % From KVB paper (division by 2 from KVB note)
        KVBmat1=repmat(KVBcrits(1,:),draws,1);
    end
    basisfncrits=finv(.95,m,dof-2+1); % Exact F for fixed-b
    basismat1=repmat(basisfncrits(1,:),draws,1);
    
    % Rejection rates
    QSrej=1-mean(FQS<QSmat1);
    NWrej=1-mean(FNW<NWmat1);
    EWPrej=1-mean(FEWP<basismat1);
    COSrej=1-mean(FCOS<basismat1);
    SSrej=1-mean(FSS<basismat1);
    KVBrej=1-mean(FKVB<KVBmat1);
    LGrej=NaN(size(SSrej));
    MCsize=[QSrej,NWrej,EWPrej,COSrej,SSrej,KVBrej];
    
    % Quantiles
    QSquantiles=quantile(FQS,.95);
    QSmat1=repmat(QSquantiles(1,:),draws,1);
    NWquantiles=quantile(FNW,.95);
    NWmat1=repmat(NWquantiles(1,:),draws,1);
    EWPquantiles=quantile(FEWP,.95);
    EWPmat1=repmat(EWPquantiles(1,:),draws,1);
    COSquantiles=quantile(FCOS,.95);
    COSmat1=repmat(COSquantiles(1,:),draws,1);
    SSquantiles=quantile(FSS,.95);
    SSmat1=repmat(SSquantiles(1,:),draws,1);
    KVBquantiles=quantile(FKVB,.95);
    KVBmat1=repmat(KVBquantiles(1),draws,1);
    
    % Compute alternatives
    delrange=0:.25:4;
    if isequal(DGP,@stochasAR)
        rho=alpha^2;
        betas=delrange*sqrt((1+rho)/(1-rho))/sqrt(T);
    elseif isequal(DGP,@loc2)
        betas=delrange/(1-alpha)/sqrt(T);
    end
    
    % Construct Storage
    QSpow=NaN(length(delrange),4);
    NWpow=NaN(length(delrange),4);
    EWPpow=NaN(length(delrange),4);
    COSpow=NaN(length(delrange),4);
    SSpow=NaN(length(delrange),4);
    KVBpow=NaN(length(delrange),1);
    FQS=zeros(draws,length(dof));
    FNW=zeros(draws,length(dof));
    FEWP=zeros(draws,length(dof));
    FCOS=zeros(draws,length(dof));
    FSS=zeros(draws,length(dof));
    FKVB=zeros(draws,1);
    
    % Test stats
    for j=1:length(delrange)
        betatest=betadev+[repmat(betas(j),draws,1),zeros(draws,1)]; % Augment estimated beta deviation from run with beta_alt to calculate power under alterntive
        for i=1:draws
            for d=1:length(dof)
                FQS(i,d)=T*betatest(i,:)*squeeze(omegastoreQS(i,d,:,:))^-1*betatest(i,:)'/m;
                FNW(i,d)=T*betatest(i,:)*squeeze(omegastoreNW(i,d,:,:))^-1*betatest(i,:)'/m;
                FEWP(i,d)=factor(d).*T*betatest(i,:)*squeeze(omegastoreEWP(i,d,:,:))^-1*betatest(i,:)'/m;
                FCOS(i,d)=factor(d).*T*betatest(i,:)*squeeze(omegastoreCOS(i,d,:,:))^-1*betatest(i,:)'/m;
                FSS(i,d)=factor(d).*T*betatest(i,:)*squeeze(omegastoreSS(i,d,:,:))^-1*betatest(i,:)'/m;
                
            end
            FKVB(i)=T*betatest(i,:)*squeeze(omegastoreKVB(i,1,:,:))^-1*betatest(i,:)'/m;
        end
        
        % Rejection Rates
        QSpow(j,:)=1-mean(FQS<QSmat1);
        NWpow(j,:)=1-mean(FNW<NWmat1);
        EWPpow(j,:)=1-mean(FEWP<EWPmat1);
        COSpow(j,:)=1-mean(FCOS<COSmat1);
        SSpow(j,:)=1-mean(FSS<SSmat1);
        KVBpow(j,:)=1-mean(FKVB<KVBmat1);
    end
    powsample=[QSpow,NWpow,EWPpow,COSpow,SSpow,KVBpow];
    
    %% Find worst case power
    grid=0:.05:4;
    samp=delrange;
    pow_gau=1-ncx2cdf(chi2inv(.95,m),m,grid.^2); % Oracle power
    powinterp=interp1(samp, powsample,grid,'linear'); % Interpolate curves
    [maxPloss]=max(repmat(pow_gau',1,length(MCsize))-powinterp); % Find maximum loss
    
end


end
