clear;

A_ej1 = [3, -0.1, -0.2; 0.1, 7, -0.3; 0.3, -0.2, 10]
b_ej1 = [7.85, -19.3, 71.4]

sol1_Seidel = metodoGaussSeidel(A_ej1, b_ej1, [0, 0, 0], 0.01)


function isDominant = diagDom(A)
    d = abs(diag(A));
    suma_fuera = sum(abs(A), 2) - d;
    isDominant = all(d > suma_fuera);
end
function isSquare = cuadrada(A)
    isSquare = size(A,1) == size(A,2);
end

function x = metodoGaussSeidel(A, b, x_inicial, epsilon)
    if ~(cuadrada(A) && diagDom(A))
        error('Debe ser un sistema compatible determinado, diagonalmente dominante.')
    end
    n = length(A(1,:))
    if ~(numel(b)==n && numel(x_inicial)==n)
        error('Las dimensiones no son correctas. ')
    end
    n = length(A(1,:))
    x = zeros(1, n);
    x_nuevo = zeros(1, n) + epsilon + 1;   %Inicializo con valor auxiliar 'lejano'
    while(all(abs(x_nuevo-x) >= epsilon))
        x = x_nuevo

        for i = 1:n
            suma_anterior = 0;
            suma_posterior = 0;
            if(i~=1)
                suma_anterior = sum(A(i,[1,i-1]).*x_nuevo([1:i-1]));
            end
            if(i~=n)
                suma_anterior = sum(A(i,[i+1,n]).*x([i+1:n]));
            end
            

            x_nuevo(i) = (b(i) - suma_anterior - suma_posterior)/A(i,i)
        end
    end
    x = x_nuevo    % Devuelve la solución
end

