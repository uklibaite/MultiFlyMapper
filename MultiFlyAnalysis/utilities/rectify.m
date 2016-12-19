function out = rectify(x)

    out = x .* heaviside(x);