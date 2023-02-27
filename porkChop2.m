function [ddates, tof, del_v] = porkChop2(rD, vD, rAi, vAi, rendezevous, mu, dateDepartMin, dateDepartMax, dateArrivalMin, dateArrivalMax, deltaVMax)
    

    %% Convert from au, au/day to km, km/sec


    rAi = rAi .*  1.495978707*10^8;
    rD  = rD  .*  1.495978707*10^8;

    vAi = vAi .* 1.495978707*10^8 ./ (24*60*60);
    vD  = vD  .* 1.495978707*10^8 ./ (24*60*60);
    
    initDate = datetime(2017, 1, 1);
    
    diff_de1 = seconds(dateDepartMax - dateDepartMin);
    diff_ar1 = seconds(dateArrivalMax - dateArrivalMin);

    timeDE1 = [seconds(dateDepartMin - initDate):diff_de1/1000:seconds(dateDepartMax - initDate)];
    timeA1  = [seconds(dateArrivalMin - initDate):diff_ar1/1000:seconds(dateArrivalMax - initDate)];
    
uop_e1_r = [];
uop_e1_v = [];
uop_s1_r = [];
uop_s1_v = [];

for (i = 1:1:length(timeDE1))
    [temp1, temp2] = universalOrbitPropogator(rD, vD, timeDE1(i));
    uop_e1_r = [uop_e1_r; temp1];
    uop_e1_v = [uop_e1_v; temp2];
end

for (i = 1:1:length(timeA1))
    [temp1, temp2] = universalOrbitPropogator(rAi, vAi, timeA1(i));
    uop_s1_r = [uop_s1_r; temp1];
    uop_s1_v = [uop_s1_v; temp2];
end

    ddates = [];
    tof = [];
    del_v = [];
    array_delv = [];
    
        for(i = 1:1:length(uop_e1_r(:, 1)))
            i
            for(j = 1:1:length(uop_s1_r))
                if (timeA1(j) - timeDE1(i)) < 0
                    continue;
                end
                  for(k = 0:1:1)  
                    [v_1, v_2, ~] = LambertCurtis(uop_e1_r(i, :), uop_s1_r(j, :), timeA1(j) - timeDE1(i), k, mu);
                    
                    if(rendezevous)
                        if (norm(v_1-uop_e1_v(i, :)) + norm(v_2-uop_s1_v(j, :))) < deltaVMax
                            
                            del_v = [del_v, norm(v_1-uop_e1_v(i, :)) + norm(v_2-uop_s1_v(j, :))];
                            tof = [tof, timeA1(j) - timeDE1(i)] ;
                            ddates = [ddates, timeDE1(i)];
                        end
                    else
                        if norm(v_1-uop_e1_v(i, :)) < deltaVMax
    
                            del_v = [del_v, norm(v_1-uop_e1_v(i, :))];
                            tof = [tof, timeA1(j) - timeDE1(i)] ;
                            ddates = [ddates, timeDE1(i)];
                            
                        end
                    end
                  end
            end

        end
   
    


tof_days   = tof/(60*60*24);
%ddt_yrs_temp    = dateDepartMin + seconds(ddates);

% ddt_yrs = [];
% 
% ddt_yrs = year(ddt_yrs_temp) + (ddt_yrs_temp - datetime(year(ddt_yrs_temp), 1, 1)) ...
% ./(datetime(year(ddt_yrs_temp), 12, 31) - datetime(year(ddt_yrs_temp), 1, 1));

% Converting b/w seconds and years made the plot ugly
%[ddq, tofq] = meshgrid(min(ddt_yrs):3.2e-4:max(ddt_yrs), min(tof_days):.11:max(tof_days));
[ddq, tofq] = meshgrid(min(ddates):10000:max(ddates), min(tof_days):.11:max(tof_days));
del_vq = griddata(ddates, tof_days, del_v, ddq, tofq);

figure
hold on
contourf(ddq, tofq, del_vq); c = colorbar; c.Label.String = 'Delta V (km/sec)';
xlabel('Departure Time from JD 2457754.5 [seconds]'); ylabel('Time of Flight [days]');


hold off



end