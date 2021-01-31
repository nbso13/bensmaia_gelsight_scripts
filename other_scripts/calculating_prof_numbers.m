%% Calculating Profilometry Numbers
indentor_diameter = 3; %mm - probably between 3mm and 5 mm
indentor_radius = indentor_diameter/2;
SA = indentor_radius^2 * pi/(10^(6)); %SA of sphere /2 , units m^2
force_per_ind = 0.08; %N/mm, from instron finger data
pressure_per_ind = force_per_ind/(1000*SA); %kPa per mm
disp(pressure_per_ind)
%how many mms to match 5kPa?
ind = 5000/pressure_per_ind;
%get
SA_mold = 0.35/(10^6); % in m sq
force = 5000*SA_mold*28;

%%

indentor_diameter = 3; %mm - probably between 3mm and 5 mm
indentor_radius = indentor_diameter/2;
SA = indentor_radius^2 * pi/(10^(6)); %SA of sphere /2 , units m^2

gel_depth = 10;%mm
depth1 = 1;
depth2 = 1.938;
force1 = 0.082; %from toastertest analysis script plot
force2 = 0.312;

E1 = ((force1/SA)/(depth1/gel_depth))/1000000 %kilo pascals
E2 = ((force2/SA)/(depth2/gel_depth))/ 1000000 %kilopascals



indentor_diameter = 5; %mm - probably between 3mm and 5 mm
indentor_radius = indentor_diameter/2;
SA = indentor_radius^2 * pi * 2/1000000; %SA of sphere /2 , units m^2
force_per_ind = 0.1; %N/mm, from instron finger data
pressure_per_ind = force_per_ind/SA; %N/m^2 per mm
SA_mold = 450/ 1000000; %m^2
pressure_range_match = [3000, 5000]; %N/m^2
ind_range = pressure_range_match./pressure_per_ind; %mm to indent to achieve pressure
ind_per_turn = 0.635; %mm per full turn of 40thread per inch screw
pressure_diff_per_turn = ind_per_turn*pressure_per_ind;
eight_turn_error = 1/8*[pressure_diff_per_turn/pressure_range_match(1), pressure_diff_per_turn/pressure_range_match(2)];

ind = 1; % mm indentation into gel
disp(strcat("pressure at indentation = ", num2str(ind*pressure_per_ind)));