function [y1,xf1,xf2] = larm_narx_1del(x1,x2,xi1,xi2)
%MYNEURALNETWORKFUNCTION neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 15-Jan-2019 19:45:23.
%
% [y1,xf1,xf2] = myNeuralNetworkFunction(x1,x2,xi1,xi2) takes these arguments:
%   x1 = 3xTS matrix, input #1
%   x2 = 2xTS matrix, input #2
%   xi1 = 3x1 matrix, initial 1 delay states for input #1.
%   xi2 = 2x1 matrix, initial 1 delay states for input #2.
% and returns:
%   y1 = 2xTS matrix, output #1
%   xf1 = 3x1 matrix, final 1 delay states for input #1.
%   xf2 = 2x1 matrix, final 1 delay states for input #2.
% where TS is the number of timesteps.

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [0;0;0];
x1_step1.gain = [2.22222222222222;2.22222222222222;2.22222222222222];
x1_step1.ymin = -1;

% Input 2
x2_step1.xoffset = [-8.82318721934812;-9.77857114320086];
x2_step1.gain = [0.113207268061184;0.107408102005448];
x2_step1.ymin = -1;

% Layer 1
b1 = [-2.5036919926575183;1.0870880345663556;-0.64741309315552609;-0.51252159739340619;-1.7861566753627984;0.48362462363560543;0.98361758484349426;3.1567455693359858;-0.032258116891000346;2.1583914392789327];
IW1_1 = [0.063568026546937009 -0.56991950438139727 -1.1582500883743085;0.014555370246841193 0.0023879908783840109 -0.047222673544136506;0.013133531540498037 -0.0134069451791833 -0.051677502455977388;-0.04810876268438722 0.055042671782482032 0.2055408793510892;-0.79074509462954767 1.657134210404736 1.3474425627936795;0.0075688194221387797 -0.01529903106410878 -0.045243164124667293;-0.070574624910736433 0.072063012877632115 0.32084999196642827;2.0200186319109839 1.4388345329280174 -0.79599067702069293;0.012939745748815259 -0.015302015751419403 -0.058851447799098532;1.1240740339198212 -1.7750436119563637 -1.4595937040837565];
IW1_2 = [2.2843968005514972 -1.9843148948671556;-0.16490228509785737 0.41386889217129902;-0.13858522866097786 0.40088536692105192;0.26076774296692945 -0.0044160445211053221;1.0313680479431659 -1.0745492909439649;-0.10491126595302593 0.33636985183603924;0.39282107964253754 0.012534953645289783;0.035452496364755087 4.8554671124083688;0.17564411784667447 -0.004127623494352162;-1.0718426136060322 1.0891402664912331];

% Layer 2
b2 = [0.09537979429275735;-0.28015175427014066];
LW2_1 = [0.0069025180098217604 0.55579360822634616 0.14315507780480735 1.0766507736347315 -0.042315618699067845 -0.41602894731209117 0.4895414523927169 0.0021839617961038013 4.0509440833330235 -0.034829583022025809;0.00060958063496309239 1.1195779124562504 1.8287473970653749 0.78001441097082913 -0.12094336696681308 1.2264544674391606 0.34936633738422185 0.0030303897633085653 0.72666412658223034 -0.11100869436710024];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = [0.113207268061184;0.107408102005448];
y1_step1.xoffset = [-8.82318721934812;-9.77857114320086];

% ===== SIMULATION ========

% Dimensions
TS = size(x1,2); % timesteps

% Input 1 Delay States
xd1 = mapminmax_apply(xi1,x1_step1);
xd1 = [xd1 zeros(3,1)];

% Input 2 Delay States
xd2 = mapminmax_apply(xi2,x2_step1);
xd2 = [xd2 zeros(2,1)];

% Allocate Outputs
y1 = zeros(2,TS);

% Time loop
for ts=1:TS
    
    % Rotating delay state position
    xdts = mod(ts+0,2)+1;
    
    % Input 1
    xd1(:,xdts) = mapminmax_apply(x1(:,ts),x1_step1);
    
    % Input 2
    xd2(:,xdts) = mapminmax_apply(x2(:,ts),x2_step1);
    
    % Layer 1
    tapdelay1 = reshape(xd1(:,mod(xdts-1-1,2)+1),3,1);
    tapdelay2 = reshape(xd2(:,mod(xdts-1-1,2)+1),2,1);
    a1 = tansig_apply(b1 + IW1_1*tapdelay1 + IW1_2*tapdelay2);
    
    % Layer 2
    a2 = b2 + LW2_1*a1;
    
    % Output 1
    y1(:,ts) = mapminmax_reverse(a2,y1_step1);
end

% Final delay states
finalxts = TS+(1: 1);
xits = finalxts(finalxts<=1);
xts = finalxts(finalxts>1)-1;
xf1 = [xi1(:,xits) x1(:,xts)];
xf2 = [xi2(:,xits) x2(:,xts)];
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
y = bsxfun(@minus,x,settings.xoffset);
y = bsxfun(@times,y,settings.gain);
y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
x = bsxfun(@minus,y,settings.ymin);
x = bsxfun(@rdivide,x,settings.gain);
x = bsxfun(@plus,x,settings.xoffset);
end
