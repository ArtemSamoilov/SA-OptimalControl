% function that sets psi on iteration i
function [dpsi] = iterpsi(t, C, psi)
    dpsi = (-C')*psi;
end