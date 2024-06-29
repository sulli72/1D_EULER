function [FPLS,FMNS] = CHARSPLIT(FLXVEC,STATE,alpha1,alpha2,alpha3)

%%% LAX FRIEDRICHS SPLITTING
% FPLS = 1/2*(FLXVEC + alpha*Q);
% FMNS = 1/2*(FLXVEC - alpha*Q);

%%% GLOBAL LAX FRIEDRICHS CHARACTERISTIC FORM OF FLUX SPLITTING
AVEC = [alpha1 alpha2 alpha3];
mag=max(abs(AVEC));
MYALPH=mag;

%%% GLOBAL APPROACH
FPLS(1,:) = 0.5*(FLXVEC(1,:) + MYALPH*STATE(1,:)); % <--- order needs to be consistent with order of eigenvectors
FPLS(2,:) = 0.5*(FLXVEC(2,:) + MYALPH*STATE(2,:));
FPLS(3,:) = 0.5*(FLXVEC(3,:) + MYALPH*STATE(3,:));

FMNS(1,:) = 0.5*(FLXVEC(1,:) - MYALPH*STATE(1,:)); % <--- order needs to be consistent with order of eigenvectors
FMNS(2,:) = 0.5*(FLXVEC(2,:) - MYALPH*STATE(2,:));
FMNS(3,:) = 0.5*(FLXVEC(3,:) - MYALPH*STATE(3,:));

%%% ROE-TYPE APPROACH -- SOMETIMES LOCAL ROE TYPE FLUX SPLITTING IS MORE
%%% ROBUST THAN GLOBAL SCALR TYPE FLUX SPLITTING

% alpha1=abs(alpha1);
% alpha2=abs(alpha2);
% alpha3=abs(alpha3);
% 
% FPLS(1,:) = 0.5*(FLXVEC(1,:) + alpha1*STATE(1,:)); % <--- order needs to be consistent with order of eigenvectors
% FPLS(2,:) = 0.5*(FLXVEC(2,:) + alpha2*STATE(2,:));
% FPLS(3,:) = 0.5*(FLXVEC(3,:) + alpha3*STATE(3,:));
% 
% FMNS(1,:) = 0.5*(FLXVEC(1,:) - alpha1*STATE(1,:)); % <--- order needs to be consistent with order of eigenvectors
% FMNS(2,:) = 0.5*(FLXVEC(2,:) - alpha2*STATE(2,:));
% FMNS(3,:) = 0.5*(FLXVEC(3,:) - alpha3*STATE(3,:));

end