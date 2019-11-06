function [var1a, var2a] =  mix_variables(fract1, fract2, fract1_new,...
    fract2_new, var1, var2)

%Check if sum of fractions is 1
if fract1 + fract2 ~= 1
    error('sum of fractions is not one')
end

if fract1_new + fract2_new ~= 1
    error('sum of new fractions is not one')
end


if fract1 > fract1_new
    var1a = var1;
    var2a = var2*fract2/fract2_new + var1*(fract2_new-fract2)/fract2_new;

    
elseif fract1 < fract1_new
    var2a = var2;
    var1a = var1*fract1/fract1_new + var2*(fract1_new-fract1)/fract1_new;
    
else %frac1 === frac1_new
    var1a = var1;
    var2a = var2;
end


