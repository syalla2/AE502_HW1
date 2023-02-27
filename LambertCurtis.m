function [v_1, v_2, orbital_elements] = LambertCurtis(r1,r2, t_f, pro, mu)
    %mu = 1.327124400189e20;

    %% Curtis Algorithm 5.2, Step 1 
    mag_r1 = norm(r1);
    mag_r2 = norm(r2);
    r1xr2 = cross(r1, r2);

    %% Curtis Algorithm 5.2, Step 2    
    del_theta = 0;
    if pro
        if(r1xr2(3) >= 0)
            del_theta = acos(dot(r1, r2)/(mag_r1 * mag_r2));
        else
            del_theta = 2*pi - acos(dot(r1, r2)/(mag_r1 * mag_r2));
        end
    else   
        if(r1xr2(3) >= 0)
            del_theta = 2*pi - acos(dot(r1, r2)/(mag_r1 * mag_r2));
        else
            del_theta = acos(dot(r1, r2)/(mag_r1 * mag_r2));
        end
    end

    %% Curtis Algorithm 5.2, Step 3
    A = sin(del_theta) * sqrt((mag_r1 * mag_r2) / (1 - cos(del_theta)));

    %% Curtis Algorithm 5.2, Step 4
    z_0 = rand(1)*.001;
    z_new = z_0;
    iter = 0;
    y_0 = 0;
    err = 10;

    while (err > 1e-08 && iter < 9)
        if(z_new > 0)
            S = ( sqrt(z_new) - sin(sqrt(z_new)) ) / ( sqrt(z_new)^3 );
            C = ( 1 - cos( sqrt(z_new) ) )/z_new;
        elseif(z_new < 0)
            S = ( sinh(sqrt(-z_new)) - (sqrt(-z_new))) / ( sqrt(-z_new)^3 );
            C = ( cosh( sqrt(-z_new) ) - 1 )/ (-z_new);
        elseif(z_new == 0)
            S = 1/6;
            C = 1/2;
        end
        
        y = mag_r1 + mag_r2 + A * (z_new*S - 1)/sqrt(C);
        if iter == 0
            y_0 = y;
        end


        F = (y/C) ^ (3/2) * S + A*sqrt(y) - sqrt(mu)*t_f;
        

        if(z_new == 0)
            F_prime = sqrt(2)/40 * y_0^(3/2) + A/8 * (sqrt(y_0) + A*sqrt(1 / (2 * y_0) ) );
        else
            F_prime = (y/C)^(3/2) * (1/(2*z_new) * (C - 1.5 * S/C) + .75 * S^2/C) + A/8 * (3*S/C*sqrt(y) + A*sqrt(C/y));
        end

        iter = iter + 1;

        z_old = z_new;
        z_new = z_old - F/F_prime;
        err = norm(z_new - z_old);
    end
    %% Curtis Algorithm 5.2, Step 5
    y = mag_r1 + mag_r2 + A * (z_new*S - 1)/sqrt(C);
    
    %% Curtis Algorithm 5.2, Step 6
    f       = 1 - y/mag_r1;
    g       = A * sqrt(y/mu);
    f_prime = sqrt(mu)/(mag_r1 * mag_r2) *sqrt(y/C) * (z_new*S - 1);
    g_prime = 1 - y/mag_r2;

    %% Curtis Algorithm 5.2, Step 7

    v_1 = 1/g * (r2 - f*r1); %km/sec
    v_2 = 1/g * (g_prime * r2 - r1); %km/sec

    %% Curtis Algorithm 5.2, Step 8
    % Curtis Algorith 4.2
    v_r = dot(r1, v_1./mag_r1);
    mag_v1 = norm(v_1);
    h = cross(r1, v_1);
    mag_h = norm(h);
    i = acos(h(3)/mag_h);
    N = cross([0, 0, 1], h);
    mag_n = norm(N);
    raan = acos(N(1)/mag_n);
    if(N(2) < 0)
        raan = 2*pi - raan;
    end
    e = 1/mu * ((mag_v1^2 - mu/mag_r1)*r1 - mag_r1*v_r*v_1);
    mag_e = norm(e);

    w = acos(dot(N/mag_n, e/mag_e));
    if(e(3) < 0)
        w = 2*pi - w;
    end
    
    theta = acos(dot(e/mag_e, r1/mag_r1));
    
    if(v_r < 0)
        theta = 2*pi - theta;
    end


    orbital_elements = [mag_h, i, mag_n, raan, mag_e, w, theta];


end