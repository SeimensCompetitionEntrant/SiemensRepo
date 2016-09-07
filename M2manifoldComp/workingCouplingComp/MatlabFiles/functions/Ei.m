function value = Ei(x)

if(x/1j > 0)
    value = cosint(x/1j) + 1j*sinint(x/1j);
elseif(x/1j < 0)
    value = cosint(x/1j) - 1j*sinint(x/1j);
else
    disp('input is 0')
    value = [];
end