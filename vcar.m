function v = vcar(d)
    global dmin dmax vmax;
    if (d <= dmin)
        v = 0;
    else
        v = vmax * log(d ./ dmin) / log(dmax / dmin) .* (d > dmin) .* (d < dmax) + vmax * (d >= dmax);
    end
end
