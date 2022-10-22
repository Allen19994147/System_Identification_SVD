%%% Hankel Matric Construction
%%% Allen Lee
function H = Hankel_Matrix(y, dim)

if((2*dim-1)>length(y))
    warning("Extend data samples or decrease matrix dimension!")
else
    H = [];
    
    for i=1:dim
        h = [y(i)];
        for j = 1:dim-1
            h = [h;y(i+j)];
        end
        H = [H h];
    end

end
