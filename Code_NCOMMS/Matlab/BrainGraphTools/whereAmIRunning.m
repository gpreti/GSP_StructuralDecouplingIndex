function location=whereAmIRunning();
%
% host machine identification for nomad users
%
% returns
%  'pleiades' if executing on any pleiades node
%  'laptop' if executing on my laptop
%  'its' if executing on the itspc31 computer
% v1.0 Dec 2004 Jonas Richiardi
% v1.0.1 2005 Aug 12 JR - fix new host name for laptop
% v1.0.2 2006 Nov 8 JR - fix hostname for its server
% v1.0.3 2008 Mar 12 JR - fix new laptop name
% v1.0.4 2008 Oct JR - ditto
% v1.0.5 2009 Sep JR - yields errors for unknown hosts
% v1.0.6 2009 Oct JR - fix for macBookPro + dynamic hostname

[returnCode, hostName]=system('hostname');

% on MacOS, hostname is not static by default. see here:
% http://excitedcuriosity.wordpress.com/2007/08/24/mac-os-x-hostname-determ
% ination/
if (strcmp(deblank(hostName),'jmac.local') || ismac) % MacBookPro laptop
    location='jmac';
elseif(strcmp(deblank(hostName),'null3')) % dell laptop
    location='laptop';
elseif(strcmp(deblank(hostName),'null')) % toshiba laptop
    location='laptop';
else
    error('Unknown host name, not sure where I am running');
end
