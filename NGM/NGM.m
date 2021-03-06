%% Description: Gradient free method
% Input:
% 1. num_iter - number of iteration of algo
% 2. t_start - start point
% 3. nu - sensitivity to stoh gradient.
% 2. mu - model constant
% 3. d - demand flow
% 4. e - accurancy
% 6. L - length of paximal path in graph
% 5. net - information about graph edges (node_o, node_d, time, flow)
% 6. mc - matrix connectivity
% 7. pairs - all posible pairs of nodes in given graph (node_o, node_d)
% Output: 
% 1. t - optimal argument
% 2. F - optimal function value

function [t, fun] = NGM(num_iter, t_start, R, nu, mu, d, L, net, mc, pairs)
    R_show = R
    for k = 1:num_iter
        k;
        u = randn(size(t_start,1),1);
        h = R/sqrt(k);
        t_start;
        fun      =  stohdm(t_start        , d, L, mu, net, mc, pairs);
        step = t_start + nu.*u;
        fun_step =  stohdm(step, d, L, mu, net, mc, pairs);
        grad = (fun_step - fun)/nu;
        t = t_start - grad*h*u;
        t = proj(net,t);
%         opt_diff = t - net(:,3);
        t_start = t;
        
    end
end

