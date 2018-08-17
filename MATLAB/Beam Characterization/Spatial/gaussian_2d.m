function [gaussian_profile, x, y] = gaussian_2d(XY, peak, waist_major, waist_minor, angle_degrees, center_x, center_y)

X = XY(:,1);
Y = XY(:,2);

a = ((cosd(angle_degrees)^2) / (waist_major^2)) + ((sind(angle_degrees)^2) / (waist_minor^2));
b = -((sind(2*angle_degrees)) / (2*waist_major^2)) + ((sind(2*angle_degrees)) / (2*waist_minor^2));
c = ((sind(angle_degrees)^2) / (waist_major^2)) + ((cosd(angle_degrees)^2) / (waist_minor^2));

gaussian_profile = peak*exp(-(a*(X - center_x).^2 + 2*b*(X - center_x).*(Y - center_y) + c*(Y - center_y).^2));



end

