function estX         = estFactor(L, Ph, Y) % estimation of factor activity
    estX              = L'/(L*L'+diag(Ph)) * Y';
end