function lambda = hazardRate_nonExtrap(x_desired, events)

% David Huberdeau, 04/28/16
%
% Compute the hazard rate function at the points/times x_desired using
% right censored flag censors.
%
% primarily for special case of skill game analysis.

global alignment_error_tolerance

if ~exist('alignment_error_tolerance', 'var')
    alignment_error_tolerance = .01;
end

% [f, x] = ecdf(events, 'censoring', censors);
[f, x] = ecdf(events);

cdf = nan(length(x_desired), 1);
for i_x = 1:length(x_desired)
    % new way... step cdf value only when it increments, but without
    % extrapolation:
    if x_desired(i_x) > max(x)
        % the last recorded falloff has occured, do NOT extrapolate
        cdf(i_x) = nan;
    else
        % have not yet reached the last falloff, so keep recording hazard
        [x_error, kmin] = min(abs(x - x_desired(i_x)));
        cdf(i_x) = f(kmin);
    end
    
    % old flawed way...
%     [x_error, kmin] = min(abs(x - x_desired(i_x)));
%     if x_error < alignment_error_tolerance
%         cdf(i_x) = f(kmin);
%     else
% %         cdf(i_x) = nan;
%         if x(kmin) == x(end) && (x(kmin) - x_desired(i_x)) < 0
%             % past end of signal
%             cdf(i_x) = f(kmin-1); % retain last value
%         elseif x(kmin) == x(end) && (x(kmin) - x_desired(i_x)) >= 0
%             % near end of signal but not past end
%             cdf(i_x) = f(kmin);
%         elseif x(kmin) ~= x(end) && kmin ~= 1 && (x(kmin) - x_desired(i_x)) < 0
%             % middle of signal; 
%             cdf(i_x) = f(kmin); %take last nearest cdf estimate
%         elseif x(kmin) ~= x(end) && kmin ~= 1 && (x(kmin) - x_desired(i_x)) >= 0
%             cdf(i_x) = f(kmin-1); %take last nearest cdf estimate
%         elseif kmin == 1
%             cdf(i_x) = f(kmin); %take estimate of first cdf
%         else
%             % not sure what is left
%             cdf(i_x) = nan;
%             warning('Should this have happened?');
%         end
%     end
end
% cdf_filt = sgolayfilt(cdf, 5, 13);
cdf_filt = cdf; % no filtering yet (sometimes filtering makes slope of cdf negative, which is impossible)
df = [0; diff(cdf_filt)./diff(x_desired)'];

cdf_filt_x1 = cdf_filt(cdf_filt < 1); 
cdf_max_x1 = max(cdf_filt_x1);
cdf_filt(cdf_filt == 1) = cdf_max_x1; % no dividing by zero

if sum(cdf_filt > 1) > 0
    error('Fatal error computing hazard rate function.')
end

% minor filtering to make H.R. functions less sharp
% lambda = sgolayfilt(df./(1 - cdf_filt), 3, 5); %window size 10% of x_des
lambda = df./(1 - cdf_filt); %this derives from the hazard rate definition
lambda(~isfinite(lambda)) = nan;
lambda(lambda < 0) = 0;


% lambda = df./(1-cdf_filt);
% lambda(~isfinite(lambda)) = nan;


