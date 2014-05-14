clc; clear;

%% Initial parameters
%Net generation parameters
len = 4;
high = 4;
save_flag = 0;
load_flag = 1;
path_to_file = '/Users/Egor/Desktop/SDModel/Alpha/g4x4.mat';

%% Generate manhattan net
%net - information about time on each arc. (origin, destination, time, flow)
[net, mc] = generate_manhatan(len, high, save_flag, load_flag, path_to_file);
pairs = making_pairs_func(len*high);

%% method part

net
L = 10;
e = 0.1;
mu = 0.1;

d = zeros(size(pairs,1),1);
d(15) = 0.2;

R = 1;

nu = 0.1;




%% test part

t_start = ones(size(net,1),1);
num_iter = 10;
t_opt = net(:,3);

%1.1 Finding R for the task

r_test_num = 10;
test_num = 10;


record_vec = zeros(r_test_num,1);
R_vec = zeros(r_test_num,1);
set_start_diff = zeros(r_test_num,1);


R_vec(1) = 1;
record = 0;

for j = 1:r_test_num
    R_vec(j)
    for i = 1:test_num
        [t, f_current] = NGM(num_iter, t_start, R_vec(j), nu, mu, d, L, net, mc, pairs);
        dist = sqrt((t - t_opt)'*(t - t_opt));
        if dist > record
            record = dist;
        end
    end
    record_vec(j) = record;
    if abs(record) > R_vec(j) && j ~= r_test_num
        R_vec(j+1) = record;
    else 
        R_vec(j+1) = R_vec(j)/2;
    end
    set_start_diff(j) = record - sqrt((t_opt - t_start)'*(t_opt-t_start));
    record = 0;
    
    
end

record_vec;

counter = 0;
for i = 1:r_test_num
    if record_vec(i) == 0
        counter = counter + 1;
    end
end

R_opt = sum(record_vec)/(r_test_num - counter)

%2. Cheking algorythm for working.

% num_iter = 500;
% 
% [t, f_opt] = NGM(num_iter, t_opt, nu, R_opt, mu, d, L, net, mc, pairs);
% opt = f_opt
% 
% while true 
%     [t, f_current] = NGM(num_iter, t_start, R_opt, nu, mu, d, L, net, mc, pairs);
%     
%     curr = f_current
%     if abs(f_current - f_opt) < e
%         break
%     else
%         num_iter = num_iter + 500
%         abs(f_current - f_opt)
%     end
% end
% num_iter_show = num_iter


%% Output part
%graph(mc, net, net(:,4), len, high)




