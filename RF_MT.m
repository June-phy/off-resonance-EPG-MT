%-------------------------------------------------------------------------
% Function:    RF_MT
% Authors:     Qianxue Shan, Weitian Chen
% Email:       qxshan@link.cuhk.edu.hk
% Version:     1.0
% Date:        2025-10-23
%
% Copyright (c) 2025 Qianxue Shan, Weitian Chen. All rights reserved.
%
% License and Usage Notice:
%   This code is provided strictly for academic and research purposes only.
%   Any commercial use, including but not limited to sale, redistribution,
%   or integration into proprietary software, is strictly prohibited without
%   explicit written permission from the authors.
%
%   Modification of this code, its header comments, or removal of this notice,
%   in whole or in part, is EXPRESSLY FORBIDDEN without prior written consent
%   from the authors.
%
%   By using, copying, or referencing this code, you agree to abide by these terms.
%   For any inquiries or requests, please contact the authors at the email above.
%
% Description:
%   Computes the saturation rate for RF-driven magnetization transfer (MT)
%   based on the specified lineshape (SuperLorentzian, Gaussian, or Lorentzian).
%
%   Inputs:
%     T2c        (s):   Transverse relaxation time constant
%     w1         (rad/s):   RF amplitude
%     dw         (rad/s):   Frequency offset(s)
%     lineshape  (str):   Lineshape type ('SuperLorentzian', 'Gaussian', or 'Lorentzian')
%
%   Outputs:
%     rfmt (Hz): Computed saturation rate for each frequency offset in `dw`
%
%   Note:
%     For 'SuperLorentzian' lineshape, the computation avoids singularities 
%     by setting a cutoff value for `dw`.
%
%-------------------------------------------------------------------------
function rfmt=RF_MT(T2c,w1,dw,lineshape)

if strcmp(lineshape,'SuperLorentzian') %SuperLorentzian
    step=1/10000;
    sup_theta = 0:step:pi/2;
    
    cutoff=10;
    for j=1:length(dw)
        if abs(dw(j)) >= cutoff  %% the superlorentz has a pole: this avoids infinities
            % see Morrsion and Henkelman 1995.  units = s.  Seems weird to
            % me.  Need to multiply by w1^2 * pi to get saturation rate.
            du=.0001;
            u=0:du:1;
            integrand2=sqrt(2/pi)*T2c./abs(3*u.^2-1) .* exp(-2*(dw(j)*T2c./abs(3*u.^2-1)).^2);
            G(j)=w1.^2.*pi.*sum(integrand2)*du;
        else
            X = [-1.1*cutoff -cutoff cutoff 1.1*cutoff];
            Y = RF_MT(T2c,w1,X,lineshape);
            G(j)=interp1(X,Y,dw(j),'spline');
        end
        rfmt=G;
        
    end

elseif strcmp(lineshape,'Gaussian')
    rfmt=w1.^2*T2c.*sqrt(pi/2).*exp(-(dw.*T2c).^2./2); %Gaussian
else
        rfmt=w1.^2*T2c./(1+(dw.*T2c).^2); %Lorentzian
end

end
% function rfmt=RF_MT(T2c,wx,dw,lineshape)
% 
% if strcmp(lineshape,'SuperLorentzian') %SuperLorentzian
%     step=1/10000;
%     sup_theta = 0:step:pi/2;
%     cutoff=10;
%     for j=1:length(dw)
%         %if abs(dw(j)) >= cutoff  %% the superlorentz has a pole: this avoids infinities
%             % see Morrsion and Henkelman 1995.  units = s.  Seems weird to
%             % me.  Need to multiply by w1^2 * pi to get saturation rate.
% %             fun = @(theta) sin(theta)*sqrt(2/pi)*T2c./abs(3*cos(theta).^2-1) .* exp(-2*(dw(j)*T2c./abs(3*cos(theta).^2-1)).^2);
% %             g = integral(fun,0,pi/2);
% %             G(j) = wx.^2.*pi.*g;
%             du=pi/400;
%             u=0:du:pi/2;
%             integrand2=sin(u).*sqrt(2/pi)*T2c./abs(3*(cos(u)).^2-1) .* exp(-2*(dw(j)*T2c./abs(3*(cos(u)).^2-1)).^2);
%             G(j)=wx.^2.*pi.*sum(integrand2)*du;
% %         else
% %             X = [-1.1*cutoff -cutoff cutoff 1.1*cutoff];
% %             Y = RF_MT(T2c,w1,X,wc,lineshape);
% %             G(j)=interp1(X,Y,dw(j)-wc,'spline');
% %         end
%         rfmt=G;
%         
%     end
% 
% elseif strcmp(lineshape,'Gaussian')
%     rfmt=wx.^2*T2c.*sqrt(pi/2).*exp(-(dw.*T2c).^2./2); %Gaussian
% else
%         rfmt=wx.^2*T2c./(1+(dw.*T2c).^2); %Lorentzian
% end
% 
% end