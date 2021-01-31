function [pressures] = ampCurve(ts_struct, pin_radius, prof_area, gel_mass, amplitudes, plot_flag)
%ampCurve takes in a ts struct and querry amplitues, and finds the
%pressures exerted by the texture on the skin for that amplitude. Area
%comes in in millimeters. pressures go out in newtons per meter sq
%(pascals). gel_mass is the mass on the gel (i.e., weight to match with
%touchsim)

forces = zeros(size(amplitudes));
for i = 1:length(amplitudes) % for every amplitude
    disp(strcat("Indenting at (mm): ", num2str(amplitudes(i))));
    new_offset = amplitudes(i) + ts_struct.offset - max(ts_struct.offset); %setting up amplitude
    new_offset(new_offset<0) = 0;
    skin_mod_plot_flag = 0;
    [~, P] = skinModel(ts_struct.shape, new_offset', pin_radius, skin_mod_plot_flag);
    total_forces = sum(sum(P));
    forces(i) = total_forces;
end

forces = forces./prof_area; %now in N/mm^2 
pressures = forces*1000000; %now in N/m^2
if plot_flag
    figure
    plot(amplitudes, pressures)
    title("Amplitude vs pressure")
    xlabel("Amplitude of indentation (mm)")
    ylabel("Pressure (N/m^2)");
    force = gel_mass*0.0098; %grams to newtons
    gel_pressure = force/(576/1000000); %576 = 32mm*18mm, area of gel surface.
    yline(gel_pressure);
end
end

