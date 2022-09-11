%----------------------------------------------------------
% Plot undeformed shape from mesh data
%
% NEED:
%   *.node : data of nodal points coordinates
%   *.element : data of element connectivity
%
% AUTHOR:
%   Ki-Tae Kim (qlsn55@gmail.com), 25, March, 2015
%----------------------------------------------------------
clear; clc;
opengl software;
prompt = 'Input the file name (without extension) \n    >';
fname = input(prompt);

fname_mesh = strcat(fname, '.node');
fname_element = strcat(fname, '.element');
unit1 = fopen(fname_mesh, 'r');
unit2 = fopen(fname_element, 'r');

% mesh data
for i = 1 : 6
    dum = fgetl(unit1);
end

tn_node = fscanf(unit1,'%d');
node = zeros(tn_node,2);

for i = 1 : 10
    dum = fgetl(unit1);
end

for i = 1 : tn_node
    n = fscanf(unit1, '%d', 1);
    node(n,1:2) = fscanf(unit1, '%f', 2);
    dumi = fscanf(unit1, '%d', 2);
end

% element data
for i = 1 : 25
    dum = fgetl(unit2);
end

tn_element = fscanf(unit2,'%d');
conn = zeros(tn_element,4);

for i = 1 : 5
    dum = fgetl(unit2);
end

for i = 1 : tn_element
    n = fscanf(unit2, '%d', 1);
    conn(n,1:4) = fscanf(unit2, '%d', 4);
%     dumi = fscanf(unit2, '%d', 6);  % for 9-nodes
    dumi = fscanf(unit2, '%d', 1);  % for 4-nodes
end

fclose(unit1);
fclose(unit2);

% plot the mesh
xdata = zeros(4,tn_element);
ydata = zeros(4,tn_element);
for i = 1 : tn_element
    for j = 1 : 4
        n = conn(i,j);
        x = node(n,1);
        y = node(n,2);
        
        xdata(j,i) = x;
        ydata(j,i) = y;
    end
end
figure; patch(xdata,ydata,'w');

% plot nodes
hold on;
for i = 1 : tn_node
    x = node(i,1);
    y = node(i,2);
    plot(x,y,'o');
end

