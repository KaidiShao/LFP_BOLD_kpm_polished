cd D:\Koopman_data\MatlabRD
startup

addpath('D:\Onedrive\ICPBR\forNikos\MatLab\DESNET');
addpath('D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys');

which E10ha1 -all
which rpgetpars -all

SES = 'E10ha1';
sesdumppar(SES,6);
sesclnadjevt(SES,6);
sesgetcln(SES,6);
sesgetblp(SES,6);
