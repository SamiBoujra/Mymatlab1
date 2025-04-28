function A_positive = eigen_decomposition(A)
    % Décompose la matrice A en valeurs et vecteurs propres
    [V, D] = eig(A);
    
    % Extraire les valeurs propres
    eigenvalues = diag(D);
    
    % Garder uniquement les valeurs propres positives
    positive_eigenvalues = eigenvalues(eigenvalues > 0);
    
    % Créer une nouvelle matrice D_positive avec les valeurs propres positives
    D_positive = diag(positive_eigenvalues);
    
    % Recomposer la matrice avec les valeurs propres positives et les vecteurs propres correspondants
    % On garde uniquement les colonnes de V qui correspondent aux valeurs propres positives
    V_positive = V(:, eigenvalues > 0);
    
    % Recomposer la matrice en utilisant les vecteurs propres et les valeurs propres positives
    A_positive = V_positive * D_positive * V_positive';
end