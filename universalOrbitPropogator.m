%% Homework 1: Problem 1
%  Write your own Universal Variable two-body orbit propagator
%  (e.g., Curtis% Algorithms 3.3 and 3.4 [Curtis, 2013])!

function [r_, v_] = universalOrbitPropogator(r0, v0, t)

%% Algorith 3.4 Curtis
mu = 1.327124400189e11;

% Find the magnitude of r0 and v0
mag_r0 = norm(r0);
mag_v0 = norm(v0);

% The radial component of velocity
v_r = dot(r0, v0)./mag_r0;

% The reciprocal alpha of the semi-major axis
alpha = 2/mag_r0 - mag_v0^2/mu;

%Initial estimate of X_0
ratio = 10;
X_0 = sqrt(mu)*norm(alpha).*t;
X_i = X_0;


%Algorithm 3.3
iter = 0;
while(norm(ratio) > 1e-08 && iter < 16)
    %Eqn 3.52 & 3.53
    z = alpha*X_i.^2;
    if(z > 0)
        S = ( sqrt(z) - sin(sqrt(z)) ) / ( sqrt(z)^3 );
        C = ( 1 - cos( sqrt(z) ) )/z;
    elseif(z < 0)
        S = ( sinh(sqrt(-z)) -  sqrt(-z) ) / ( sqrt(-z)^3 );
        C = ( cosh( sqrt(-z) ) - 1 )/ (-z);
    elseif(z == 0)
        S = 1/6;
        C = 1/2;
    end
    
    F_X = (mag_r0 .* v_r)./sqrt(mu) .* X_i.^2 .* C + (1 -alpha*mag_r0).*X_i.^3.*S + mag_r0.*X_i - sqrt(mu).*t;
    Fdot_X = (mag_r0 .* v_r)./sqrt(mu) .*X_i.*(1 - alpha*X_i^2*S) + (1-alpha*mag_r0).*X_i^2.*C + mag_r0;

    ratio = F_X./Fdot_X;
    X_i = X_i - ratio;
    iter = iter + 1;
end
X = X_i;

f = 1 - X.^2/mag_r0 * C;
g = t - 1/sqrt(mu) * X.^3 * S;

r_ = f*r0 + g*v0;

f_dot = sqrt(mu) / (mag_r0*norm(r_)) * (alpha*X.^3*S - X);
g_dot = 1 - X.^2/norm(r_) * C;

v_ = f_dot*r0 + g_dot*v0;
end