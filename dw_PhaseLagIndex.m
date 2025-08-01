% Given a multivariate data, returns directed weighted phase lag index matrix
% dwPLI(ch1, ch2) : 
% if it is greater than 0, ch1->ch2
% if it is less than 0, ch2->ch1

% 10.28.2022. modified by Joon-Young Moon to improve speed


function dwPLI=dw_PhaseLagIndex(X)

ch=size(X,2); % Column should be channel
dwPLI=zeros(ch,ch);

%%%%%% Hilbert transform and computation of phases
a_sig=hilbert(X); % Analytic signal

%phi1=angle(a_sig); % Instantaneous phase
%mag1=abs(a_sig); % Instantaneous magnitude

for ch1=1:ch-1
    for ch2=ch1+1:ch
        
        %%%%%% Weighted phase lag index
       
%         % Algorithm #1
%         PDiff=phi1(:,ch1)-phi1(:,ch2); % phase difference
%         signed=sign(sin(PDiff)); % sign of phase difference
%         imagcomp=imag(a_sig(:,ch1).*conj(a_sig(:,ch2))); % imaginary component
%         imagmag=abs(imagcomp); % imaginary component's magnitude      
%         nom=mean(signed.*imagmag);
%         denom=mean(imagmag);                
%         dwPLI(ch1,ch2)=nom/denom;

%         % Algorithm #2
%         PDiff=mag1(:,ch1).*exp(1i.*phi1(:,ch1)).*mag1(:,ch2).*exp(-1i.*phi1(:,ch2));
%         dwPLI(ch1,ch2)=mean(abs(PDiff).*sin(angle(PDiff)))./mean(abs( abs(PDiff).*sin(angle(PDiff))));

        % Alogorithm #3: fastest
        c_sig=a_sig(:,ch1).*conj(a_sig(:,ch2));
        dwPLI(ch1,ch2)=mean(imag(c_sig))./mean(abs(imag(c_sig)));

        dwPLI(ch2,ch1)=-dwPLI(ch1,ch2);
    end
end

% Alternative algorithm #4: speed is about the same

% for ch1=1:ch
%     c_sig=repmat(a_sig(:,ch1),1,ch-ch1).*conj(a_sig(:,ch1+1:end)); % conjugate product to prevent phase jump
%     dwPLI(ch1,ch1+1:end)=mean(imag(c_sig))./mean(abs(imag(c_sig)));
%     dwPLI(ch1+1:end,ch1)=(-dwPLI(ch1,ch1+1:end))';
% end



% figure;plot(100*signed);hold on;
% plot(imagcomp,'r');
% plot(X(:,1),'g')