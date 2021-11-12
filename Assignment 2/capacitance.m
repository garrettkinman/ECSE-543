function Capacitance = capacitance(filename1, filename2, filename3)
    clc;
    
    % zero initialize
    W = 0;
    U = zeros(1, 4);

    % potential between inner and outer
    V = 15 - 0;

    % use potentials from SIMPLE2D_M and S from Question 1
    potentials = SIMPLE2D_M(filename1, filename2, filename3);
    S = [1, -0.5, 0, -.5; -0.5, 1, -0.5, 0; 0, -0.5, 1, -0.5; -0.5, 0, -0.5, 1];
    
    % e0 = epsilon naught
    e0 = 8.854e-12;
    
    for i = 1:length(potentials)
        if potentials(i, 2) < 0.1 && potentials(i, 3) < 0.1
            U(1) = potentials(i, 4);
            U(2) = potentials(i + 1, 4);
            U(3) = potentials(i + 7, 4);
            U(4) = potentials(i + 6, 4);
            W = W + (0.5 * U * S * transpose(U));
        end
    end
    Capacitance = 8 * W * e0 / V^2;
    return
end