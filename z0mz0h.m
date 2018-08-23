function [z0h]=z0mz0h(zh,nu,ustr,tstr)
%    ##################################################################
%   ##################################################################
%   ######                                                      ######
%   ######                    SUBROUTINE z0mz0h                 ######
%   ######                                                      ######
%   ######                                                      ######
%   ##################################################################
%   ##################################################################

%   input
%      real zh        ! reference level of air temperature (m)
%      real ustr      ! frictional velocity
%      real tstr      ! =-H/(rhoair*cp*ustr)
%      real nu        ! kinematic viscousity 
%   output
%      real z0h	   ! thermal roughness length

      z0h = 70 * nu ./ ustr .* exp(-7.2*sqrt(ustr).*sqrt(sqrt(abs(-tstr))));  %Yang et al. 2004
      z0h = min(zh/10,max(z0h,1.0E-10));