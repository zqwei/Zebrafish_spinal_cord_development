% 
% F test for significant power across frequency
% 
% Null hypothesis: the power spectrum is white.
% 
% Original code from Chronux -- ftestc
% 
% Input:
% sig         -- signal (1 x nTime points)
% fs          -- sampling rate
% NW          -- time-bandwidth product for tapers
% K           -- # of independent tapers to average over (< 2*NW), matlab
%                by default using K = 2*NW - 1
% pval        -- pvalue for F test
% fwindow     -- searching frequency range default: (0Hz to Nyquist freq: fs/2)

% dT = window length
% dS = window step size
% pval = criterion for F-test
%
%
% Documented and further developed by:
% Ziqiang Wei
% version 1.0
%

function SSF       = MTMFTest(sig, fs, NW, K, dT, dS, pval, fwindow)

    if nargin < 8
        fwindow    = [0 fs/2];
    end

%     dT2            = round(dT*fs);
%     dS2            = round(dS*fs);

    Nsig           = length(sig);

    [tapers, ~]    = dpsschk([NW K], Nsig, fs); % generate dpss sequence

    params         = [];
    params.tapers  = tapers;
    params.Fs      = fs;
    params.pad     = 0;
    params.fpass   = fwindow;


%     kk=ceil((length(d)-dT2+1)/dS2);

%     pos            = 1:dS2:dS2*kk;
    [f,~]      = getfgrid(params.Fs, max(2^(nextpow2(Nsig)+params.pad), Nsig), params.fpass);
    
    [Fval,A,f,sig,sd] = ftestc(data,params,pval,plt);
    
%     dim1 = length(f);
%     Fval=zeros(dim1,1, 'single');
%     
%     
%     A=zeros(dim1,kk, 'single');
%     [~,~,f,sig,~] = ftestc(d(1:(1+dT2-1)),params,pval/dT2,'y');
%     
%     events_cell = cell(size(A,2),1);
%     
%     t=(0:(size(A,2)-1))*dS2/fs;
%     
% %     parfor k=1:kk
% %         [Fval,A(:,k),~,~,~] = ftestc(d(pos(k):(pos(k)+dT2-1)),params,pval/dT2,'y');
% %         disp(Fval)
% %         fmax=crx_findpeaks(Fval,sig); %this function name is a hack. chronux 'findpeaks' conflicts with Matlab 'findpeaks'.
% %         %I have renamed the chronux function as crx_findpeaks and changed this line too.
% %         %This means this code is incompatible with the public version of chronux.
% %         %Users must use our version. Future versions of chronux are expected to
% %         %fix this namespace conflict, which will require rewrite of this line.
% %         events_cell{k}=[repmat(t(k)+dT/2,length(fmax(1).loc),1) f(fmax(1).loc)'];
% %     end
%     t = t +dT/2;
% 
%     events = cell2mat(events_cell);
%     events(:,1) = round(events(:,1)*fs);
% 
% 
%     %SSF.d=d;
%     %SSF.fs=fs;
%     %SSF.NW=NW;
%     %SSF.K=K;
%     %SSF.dT=dT;
%     %SSF.dS=dS;
%     %SSF.pval=pval;
%     SSF.t=round(t.*fs);%return time in sample units
%     SSF.f=f;
%     SSF.A=A;
%     SSF.events=events;
