% input X : data (time x channel)
% output dPLI: dPLI
% Given a multivariate data, returns phase lag index matrix
% PLI(ch1, ch2) : 
% if it is greater than 0, ch1->ch2
% if it is less than 0, ch2->ch1


function dPLI=d_PhaseLagIndex2(X)

ch=size(X,2); % column should be channel
dPLI=zeros(ch,ch);

%%%%%% Hilbert transform and computation of phases
% for i=1:ch
%     x=X(:,i);
%     %     phi0=angle(hilbert(x));  % only the phase component
%     %     phi1(:,i)=unwrap(phi0);  % smoothing
%     phi1(:,i)=angle(hilbert(x));
% end
H=hilbert(X);

%phi1=angle(H);


for ch1=1:ch-1
    for ch2=ch1+1:ch
        %%%%%% phase lage index
      %  PDiff=phi1(:,ch1)-phi1(:,ch2); % phase difference
      %  dPLI(ch1,ch2)=nanmean(sign(sin(PDiff))); % only count the sign of the asymmetry (whether the phase difference is plus or minus)
        
        HDiff=H(:,ch1).*conj(H(:,ch2)); % conjugate product : alternative method, faster speed
        dPLI(ch1,ch2)=nanmean(sign(imag(HDiff))); % only count the sign of the asymmetry (whether the phase difference is plus or minus)
        
        dPLI(ch2,ch1)=-dPLI(ch1,ch2);
    end
end




