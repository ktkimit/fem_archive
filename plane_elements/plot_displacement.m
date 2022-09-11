%----------------------------------------------------------
% Plot displements computed
%
% NEED:
%   *.ugridx : x points
%   *.ugridy : y points
%   *.ux : horizontal displacements computed
%   *.uy : vertical displacements computed
%
% AUTHOR:
%   Ki-Tae Kim (qlsn55@gmail.com), 25, March, 2015
%----------------------------------------------------------
clear; clc;
opengl software;
prompt = 'Input the file name (without extension) \n    >';
fname = input(prompt);

fname_x = strcat(fname, '.ugridx');
fname_y = strcat(fname, '.ugridy');
fname_ux = strcat(fname, '.ux');
fname_uy = strcat(fname, '.uy');
unit1 = fopen(fname_x, 'r');
unit2 = fopen(fname_y, 'r');

tn_element = fscanf(unit1, '%d', 1);
nx = fscanf(unit1, '%d', 1);
ny = fscanf(unit1, '%d', 1);

dum = fgetl(unit2);

leng1 = tn_element*ny;
leng2 = nx;

X = zeros(leng1,leng2);
Y = zeros(leng1,leng2);

for i = 1 : leng1
    X(i,:) = fscanf(unit1, '%f', leng2);
    Y(i,:) = fscanf(unit2, '%f', leng2);
end

UX = load(fname_ux);
UY = load(fname_uy);

fclose(unit1);
fclose(unit2);

fig1 = figure; hold on;
fig2 = figure; hold on;
for i = 1 : tn_element
    st = (i - 1)*ny + 1;
    en = i*ny;
    
    figure(fig1);
    mesh(X(st:en,:), Y(st:en,:), UX(st:en,:));
    figure(fig2);
    mesh(X(st:en,:), Y(st:en,:), UY(st:en,:));
end
