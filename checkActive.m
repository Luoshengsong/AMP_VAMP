function No_MD =  checkActive(xhat, Active_Set, Q)
    % Q: the cardinality of container-set for those active users ...
    % e.g. there exist 4 active users, we derive a container-set whose cardinality is 8(Q) ...
    % and then, see whether those 4 active users locates on the container-set.
    xhat_abs = abs(xhat);
    [~,SortIndice] = sort(xhat_abs ,'descend');
    container_set = SortIndice(1:Q);
    No_MD = length(Active_Set) - length( intersect(Active_Set, container_set) );
end