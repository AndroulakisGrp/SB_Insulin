%% This function is used to find the period of oscillatory data by 
... fitting a cosine curve to the data using a non-linear regression model

% The regression model (a number of cosines can be used for more complex
... waveforms)
    
% Function inputs: Initial guesses for
... Mesor,Amplitude,Period,Acrophase,t,y,options in that order
    ... t is the regressor varible
    ... y is the dependent variable
    ... varargin=[options,disp]
    ...options = 1 (default) robust regression
        ... options = 0 assumes normally distributed errors
    ... if disp = 1 Display plot, else do not.
    ... Note: Even when model fits correctly,
    ...acrophase result might symmertrically out of phase
    
% Function Ouptuts
    ... beta are the fitted parameter values
    ... residuals between the fitted model and input data
    ... J is the Jacobian from the regression analysis, which can be used
        ... for further hypothesis testing on parameters
    ... Cov is the covariance matrix
    
    

function [beta,residuals,J,Cov]=Single_Cosinor(Mesor,Amplitude,Period,Acrophase,t,y,varargin)

p0=[Mesor,Amplitude,Period,Acrophase]; % Initial guess vector

% Sinusoid model
sine_model=@(p,t)(p(1)+p(2)*cos(2*pi*(t)/p(3)-2*pi*p(4)/p(3)));

if nargin>8
    fprintf('%d',nargin)
    error('The function accepts only 6 or 7 inputs.You have  Open function for more details')
end

switch nargin
    case 8
        disp=varargin{1};
        options=varargin{2};
        
    case 7
        disp=varargin;
        options=1;
        
    case 6
        disp=1;
        options=1;
        
        
end

rng('default') % for reproducibility

if options==1
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
else
    opts=[];
end

% Non-linear regression
[beta,residuals,J,Cov]=nlinfit(t,y,sine_model,p0,opts);
if disp==1
    plot(t,y,t,sine_model(beta,t))
    legend('Input Data','Fitted Model')
end
end

