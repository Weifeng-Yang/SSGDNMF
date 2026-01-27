%% Calculate the value of the objective function
function loss=compute(var,ngmar)
    nga=var{1};
    for i=2:length(var)
        nga=nga*var{i};
    end
    loss=0.5*norm(ngmar-nga,'fro')^2;
end