%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [yes, cancelled, raIfRemembered]=askYesOrNo(msg, title, where, defaultIsYes, rememberId)
if nargin<4
    defaultIsYes=true;
    if nargin<2
        title=[];
    end
end
raIfRemembered=[];
if isempty(title)
    title= 'Please confirm...';
end
if nargin>2
    if ~isstruct(msg)
        m.msg=msg;
        m.where=where;
        if nargin>4
            m.remember=rememberId;
        end
        msg=m;
    else
        msg.where=where;
    end
end
if defaultIsYes
    dflt='Yes';
else
    dflt='No';
end
if nargout>1
    [~,yes,cancelled, raIfRemembered]=questDlg(msg, title, 'Yes', 'No', ...
        'Cancel', dflt);
else
    [~,yes,cancelled, raIfRemembered]=questDlg(msg, title, 'Yes',...
        'No', dflt);
end
end