function [FT,w] = cGrid(traj,sz)

    tr = traj;

    for int = 1:size(tr,1)
        for int2=1:size(tr,2)
            tr(int,int2) = tr(int,int2) + 0.00001*randn(1,1) + 0.00001*randn(1,1)*1i;
        end
    end

    w = voronoidens(tr(:));

    w = reshape(w,[size(tr,1) size(tr,2)]);

    N = [sz,sz];
    [xx,yy] = meshgrid(linspace(-1,1,N(1)));
    ph = double(sqrt(xx.^2 + yy.^2)<1);
    FT = NUFFT(traj,w,ph, 0,N, 2);

    
end

