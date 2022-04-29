function cnmfe_yc_deleteNeurons_locranked_v3(obj, ind, C2, Cn,pnr, folder_nm, I)
%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.
%   folder_nm: string, the folder to output images neuron by neuron.
%   I: the maximum z-projection of the video

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
% obj=neuron;
% ind=[1:size(obj.C_raw,1)];
% C2= obj.C_raw;

% TA 4/25/19 Replace empty cells with numbers 
for i = 1:ind
    if isempty(obj.Coor{i})
        obj.Coor{i} = [1 1; 2 2]; 
    end
end

if ~exist('ind', 'var') || isempty(ind)
    % display all neurons if ind is not specified.
    ind = 1:size(obj.A, 2);
elseif ind==-1
    ind = size(obj.A,2):-1:1;
end
if ~exist('C2', 'var'); C2=[]; end

if exist('folder_nm', 'var')&&(~isempty(folder_nm))
    % create a folder to save images
    save_img = true;
    cur_cd = cd();
    if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
    else
        fprintf('The folder has been created and old results will be overwritten. \n');
    end
    cd(folder_nm);
else
    save_img = false;
end

% obj.delete(sum(obj.A>0, 1)<max(obj.options.min_pixel, 1));

Amask = (obj.A>0);
ind_trim = false(size(ind));    % indicator of trimming neurons
ind_del = false(size(ind));     % indicator of deleting neurons
ctr = obj.estCenter();      %neuron's center

gSiz = obj.options.gSiz;        % maximum size of a neuron

% time
T = size(obj.C, 2);
t = 1:T;
if ~isnan(obj.Fs)
    t = t/obj.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

allA=obj.A;
for i=1:size(obj.A,2);
    allA(:,i)=allA(:,i)./max(allA(:,i));
end
meanA=max(allA,[],2);
% keyboard();
%% rank centroids based on spatial location
centroidXYsum=[];
for i=1:length(obj.Coor)
    try
    [ geom] = polygeom( obj.Coor{i}(1,:),obj.Coor{i}(2,:) ) ;
    
    end
    centroidXYsum=[centroidXYsum (geom(2)*2+geom(3)*2)*0.5];
%     keyboard();   
%     centroid{i}=geom(2:3);
end


[~,irerank]=sort(centroidXYsum);
% keyboard();
%% start viewing neurons
figure('position', [100, 100, 1024, 512]);
m=1;
while and(m>=1, m<=length(ind))
    %% full-frame view
    n=irerank(m);
    subplot(231);
%         obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
    obj.image(meanA);
%     colormap winter
    hold on
    for m2=1:length(ind);
        n2=irerank(m2);
        cont = medfilt1(obj.Coor{n2}')';
%         keyboard();
        if ind_del(m2)==0;
            if m2<m;
                plot(cont(1,:),cont(2,:),'r');
            else if m2==m;
                    
                    plot(cont(1,:),cont(2,:),'k');
                else
                    
                    plot(cont(1,:),cont(2,:),'y');
                end
            end
        end
    end
    hold off
clear    axis equal; axis off;
    if ind_del(m)
        title(sprintf('Neuron %d', ind(m)), 'color', 'r');
    else
        title(sprintf('Neuron %d', ind(m)));
    end
    %% plot another figure showing the golabal distribution
    subplot(232);
%     try
        imagesc(Cn.*pnr);
%         imagesc(pnr);
%         colormap winter
%         keyboard();
%     catch
%         obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
%     end
%     imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
    hold on
    cont = medfilt1(obj.Coor{n}')';
    plot(cont(1,:),cont(2,:),'r');
    hold off
    
    %% zoomed-in view
    subplot(233);
    obj.image(obj.A(:, ind(n)).*Amask(:, ind(n))); %
    %     imagesc(reshape(obj.A(:, ind(m)).*Amask(:,ind(m))), obj.options.d1, obj.options.d2));
    axis equal; 
    axis off;
    x0 = ctr(ind(n), 2);
    y0 = ctr(ind(n), 1);
    try
    xlim(x0+[-gSiz, gSiz]*2);
    ylim(y0+[-gSiz, gSiz]*2);
    catch
        keyboard();
    end
    
    %% temporal components
    subplot(2,3,4:5);cla;
    if ~isempty(C2)
        plot(t, C2(ind(n), :)*max(obj.A(:, ind(n))), 'linewidth', 2); hold on;
        plot(t, obj.C(ind(n), :)*max(obj.A(:, ind(n))), 'r');
        plot(t, obj.C_raw(ind(n), :)*max(obj.A(:, ind(n))), 'b');
    else
        
        plot(t, obj.C(ind(n), :)*max(obj.A(:, ind(n))));
    end
    xlim([t(1), t(end)]);
    xlabel(str_xlabel);
    
    %% TA addition 12/12/19 to overlay selected neuron on maximum z-projection of the actual field of view.
%     subplot(3,2,5); imshow(I); hold on;
%     for i = 1:size(obj.Coor,1)
%         plot(obj.Coor{i}(1,:), obj.Coor{i}(2,:), 'y-');
%         hold on;
%         plot(cont(1,:),cont(2,:),'r');
%     end
%     hold off;
    
    subplot(2,3,6); imshow(I); hold on;
    plot(cont(1,:),cont(2,:),'r');
    hold off;
    
    %% save images
    if save_img
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
        m = m+1;
    else
        fprintf('Neuron %d, keep(k, default)/delete(d)/split(s)/trim(t)/trim cancel(tc)/delete all(da)/backward(b)/end(e):    ', ind(m));
        
        temp = input('', 's');
        if temp=='d'
            ind_del(m) = true;
            m = m+1;
        elseif strcmpi(temp, 'b')
            m = m-1;
        elseif strcmpi(temp, 'da')
            ind_del(m:end) = true;
            break;
        elseif strcmpi(temp, 'k')
            ind_del(m) = false;
            m= m+1;
        elseif strcmpi(temp, 's')
            try
                subplot(333);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                tmpA = obj.A(:, ind(n));
                obj.A(:, end+1) = tmpA.*tmp_ind(:);
                obj.C(end+1, :) = obj.C(ind(n), :);
                obj.A(:, ind(n)) = tmpA.*(1-tmp_ind(:));
                obj.S(end+1, :) = obj.S(ind(n), :);
                obj.C_raw(end+1, :) = obj.C_raw(ind(n), :);
                obj.P.kernel_pars(end+1, :) = obj.P.kernel_pars(ind(n), :);
            catch
                fprintf('the neuron was not split\n');
            end
        elseif strcmpi(temp, 't')
            try
                subplot(333);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                Amask(:, ind(n)) = tmp_ind(:);
                ind_trim(m) = true;
            catch
                fprintf('the neuron was not trimmed\n');
            end
        elseif strcmpi(temp, 'tc')
            Amask(:, ind(n)) = (obj.A(:, ind(n)) > 0);
            ind_trim(m) = false;
        elseif strcmpi(temp, 'e')
            break;
        elseif ~isnan(str2double(temp))
            m = m + floor(str2double(temp));
            fprintf('jump to neuron %d / %d\n', m, length(ind));
        else
            m = m+1;
        end
    end
end
if save_img
    cd(cur_cd);
else
    obj.A(:, irerank(ind_trim)) = obj.A(:,irerank(ind_trim)).*Amask(:, irerank(ind_trim));
    obj.delete(irerank(ind_del));
    %     obj.Coor = obj.get_contours(0.9);plot(neuron.Coor{2})
end
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
obj.Coor = obj.get_contours(0.8); % energy  within the contour is 80% of the total 


end

