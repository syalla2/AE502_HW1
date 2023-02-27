    %% Homework 1 Script
    mu = 1.327124400189e11;

    %%Earth
    rE = [-1.796136509111975 * 10^-1, 9.667949206859814 * 10^-1, -3.668681017942158 * 10^-5];
    vE = [-1.720038360888334 * 10^-2, -3.211186197806460 * 10^-3, 7.927736735960840 * 10^-7];

    %%Destination
    %1I/’Oumouamoua
    r1i = [3.515868886595499 * 10^-2, -3.162046390773074, 4.493983111703389];
    v1i = [-2.317577766980901 * 10^-3, 9.843360903693031 * 10^-3, -1.541856855538041 * 10^-2];
    
    %2I/Borisov
    r2i = [7.249472033259724, 14.61063037906177, 14.24274452216359];
    v2i = [-8.241709369476881 * 10^-3, -1.156219024581502 * 10^-2, -1.317135977481448 * 10^-2];

    dateDepartMin1   = datetime(2017, 1, 1);
    dateDepartMax1   = datetime(2018, 1, 1);

    dateArrivalMin1  = datetime(2017, 9, 1);
    dateArrivalMax1  = datetime(2019, 2, 1);


    dateDepartMin2   = datetime(2017, 1, 1);
    dateDepartMax2   = datetime(2020, 8, 1);

    dateArrivalMin2  = datetime(2019, 6, 1);
    dateArrivalMax2  = datetime(2022, 2, 1);


    porkChop2(rE, vE, r1i, v1i, 0, mu, dateDepartMin1, dateDepartMax1, dateArrivalMin1, dateArrivalMax1, 20);
    porkChop2(rE, vE, r1i, v1i, 1, mu, dateDepartMin1, dateDepartMax1, dateArrivalMin1, dateArrivalMax1, 50);
    porkChop2(rE, vE, r2i, v2i, 0, mu, dateDepartMin2, dateDepartMax2, dateArrivalMin2, dateArrivalMax2, 20);
    porkChop2(rE, vE, r2i, v2i, 1, mu, dateDepartMin2, dateDepartMax2, dateArrivalMin2, dateArrivalMax2, 60);
