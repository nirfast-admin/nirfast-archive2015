function [J] = build_jacobian_cw_fl(mesh,data,omega)

% J = build_jacobian_cw_fl(mesh,data,omega)
%
% Used by jacobian_fl, builds the jacobian matrix
%
% mesh is the input mesh (variable)
% data is the data
% omega is frequency

[ncol,junk] = size(mesh.nodes);
[nrow] = length(find(mesh.link(:,3)~=0));
[nsd, msd] = size(mesh.link);

J.complexm = zeros(nrow,2*ncol);

% define parameters gamma and tau
if isfield(mesh,'gamma') == 0
    mesh.gamma = (mesh.eta.*mesh.muaf)./(1+(omega.*mesh.tau).^2);
end
f_gamma = complex(1,-omega.*mesh.tau);
f_tau = complex(0,-omega.*mesh.gamma);

k = 1;
for i = 1 : nsd
    if mesh.link(i,3) == 1
        if mesh.dimension == 2
            
            aphim_n=conj(data.aphim(:,jj));
            dphix_n=data.phix(:,i);
            
            % Calculate the gamma part here
            J.complexm(k,1:end/2) = ...
                IntFG(mesh.nodes(:,1:2),...
                sort(mesh.elements')',...
                mesh.element_area,...
                conj(data.aphim(:,mesh.link(i,2))),...
                (data.phix(:,mesh.link(i,1))).*f_gamma);
            
            % Extract log amplitude
            J.completem(k,:) = ...
                real(J.complexm(k,1:end/2)./data.complexm(k));
            
        elseif mesh.dimension == 3
            
            % Calculate the gamma part here
            J.complexm(k,1:end/2) = ...
                intfg_tet4(mesh.nodes(:,1:2),...
                sort(mesh.elements')',...
                mesh.element_area,...
                conj(data.aphim(:,mesh.link(i,2))),...
                (data.phix(:,mesh.link(i,1))).*f_gamma);
            
            %  Extract log amplitude
            J.completem(k,:) = ...
                real(J.complexm(k,1:end/2)./data.complexm(k));
            
        end
        k = k + 1;
    end
end
