function V=simplex(N)
% create regular simplex with N vertecies and N-1 dimentional space
% V - vertecies coordinates
% size(V)=[N-1 N]
% V(dimention_number,vertex_number)
% algorithm: http://en.wikipedia.org/wiki/Simplex#Cartesian_coordinates_for_regular_n-dimensional_simplex_in_Rn
% The coordinates of the vertices of a regular n-dimensional simplex can be obtained from these two properties,
%1)For a regular simplex, the distances of its vertices to its center are equal.
%2) The angle subtended by any two vertices of an n-dimensional simplex
% through its center is arccos(-1/(N-1))
V=zeros(N-1,N);


for nc=1:N-1
    V(nc,nc)=sqrt(1-sum([V(1:nc-1,nc); V(nc+1:end,nc)].^2)); % because |V(nc,:)|=1 (sum(V(nc,:).^2)=1)
    % here all V(nc,:) known, it has 1:nc - non-zeros and nc+1:end -zeros
    % now need to find V(nc,nc+1),V(nc,nc+2),...,V(nc,end) unding arccos(-1/(N-1)):
    % u(1)*v(1)+u(2)*v(2)+...u(nc-1)*v(nc-1)+u(nc)*v(nc)=-1/(N-1)
    % v(nc)=-(u(1)*v(1)+u(2)*v(2)+...+u(nc-1)*v(nc-1)+1/N)/u(nc)
    % u is V(:,nc)
    % v can be any of V(:,nc) V(:,nc+1) ... V(:,end)
    for nc1=nc+1:N
        V(nc,nc1)=-(sum(V(1:nc-1,nc).*V(1:nc-1,nc1))+1/(N-1))/V(nc,nc);
    end
end

% to check:
% see V'*V - symmetrical with 1 at diagonal and -1/(N-1) in rest
% mean(V,2) vector close to zero vector
