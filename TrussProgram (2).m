%% HOUSEKEEPING
close all;
clear all;

%% GET INPUT FILE
[file, path] = uigetfile('*.mat', 'Select the input file');
fname = strcat(path, file);
load(fname);

%% DEFINE UNITS
units = 'N'; % units of the load in matrix L, change to match problem before running code

%% CALC MATRIX A
c = size(C,2);
r = size(C,1)*2;
A = zeros(r, c); % create matrix A

S = [Sx; Sy];
A = [A,S];

m_length = zeros(c, 1); % create matrix for lengths of members

for m = 1:c
    count = 1;
    for j = 1:(r/2)
        jy = j+(r/2);
        if C(j,m) == 1
            tempC = C(:,m);
            I = find(tempC);
            d1 = [X(I(1)), Y(I(1))];
            d2 = [X(I(2)), Y(I(2))];

            if count == 1
                d = d2 - d1; 
                rn = norm(d);
                m_length(m,1) = rn; % find length of member to calc Pcrit later
                count = count + 1;
            elseif count == 2
                d = d1 - d2;
                rn = norm(d);
            end

            A(j,m) = d(1)./rn; % add to A
            A(jy,m) = d(2)./rn;
        end
    end
end

%% CALC MATRIX T
Ainv = inv(A);
T = Ainv*L;

%% LIVE LOAD ANALYSIS: FINDING CRITICAL MEMBER(S) AND MAXIMUM THEORETICAL LOAD
Wl = L(find(L));
for tm = 1:size(T,1) - 3
    R(tm,1) = T(tm)/Wl;
end

ml2 = m_length.^2; % calculate Pcrit from m_length for each member
mlm1 = ml2.^-1;
Pcrit = 2945*(mlm1);
U_Pcrit = 54./m_length; % calculate uncertainty of Pcrit from m_length for each member

wfail = -Pcrit./R; % calculate weight of failure for each member
U_wfail = -U_Pcrit./R; % calculate uncertainty of weight of failure for each member
wfail_min = min(wfail(wfail > 0)); % max theoretical load of critical member(s)

count2 = 1;
for w = 1:size(wfail) % display critical member(s)
    if wfail(w) == wfail_min
        if count2 == 1 
            U_wfail_min = U_wfail(w);
            count2 = count2 + 1;
        end
        crit_mem_out = strcat('Critical member:', 32, num2str(w)); 
        disp(crit_mem_out) % suppress this to get proper output display for the report
        
    end
end

wfail_out = strcat('Maximum theoretical load:', 32, num2str(abs(wfail_min)), 32, units, 32, '±', 32, num2str(abs(U_wfail_min)), 32, units);
disp(wfail_out) % suppress this to get proper output display for the report
%% COST CALCULATIONS
total_length = sum(m_length);
total_joints = r/2;
C1 = 10; % 10 dollars per joint
C2 = 1; % 1 dollar per inch
cost = C1*(total_joints) + C2*(total_length);

ratio = wfail_min/cost;

%% OUTPUT
tdy = datetime('today');
tdy.Format = 'MM/dd/yyyy';
tdy_char = char(tdy);
head_out = strcat('/% EK301, Section A3, Group CoPaZo: Perkins Z., McDowell C., Koutsoukos P.,', 32, tdy_char);
disp(head_out)

load_out = strcat('Load:', 32, num2str(Wl), 32, units);
disp(load_out)

mem_for_out = strcat('Member forces in', 32, units, ':');
disp(mem_for_out)
for m = 1:(size(T, 1) - 3)
    if T(m) > 0 
        toc = '(T)';
    elseif T(m) < 0
        toc = '(C)';
        T(m) = -T(m);
    elseif T(m) == 0
        toc = [];
    end
    f_out = strcat('m', num2str(m), ':', 32, num2str(T(m)), 32, toc);
    disp(f_out)
end

rea_for_out = strcat('Reaction forces in', 32, units, ':');
disp(rea_for_out)
sx1_out = strcat('Sx1:', 32, num2str(T(size(T,1)-2)));
sy1_out = strcat('Sy1:', 32, num2str(T(size(T,1)-1)));
sy2_out = strcat('Sy2:', 32, num2str(T(size(T,1))));
disp(sx1_out)
disp(sy1_out)
disp(sy2_out)

cost_out = strcat('Cost of truss:', 32, '$',num2str(cost));
disp(cost_out)

ratio_out = strcat('Theoretical max load/cost ratio in', 32, units, '/$:', 32, num2str(ratio));
disp(ratio_out)
