%Giulia 2019 
%% create randomized graph preserving degree
function Aran=rand_graph(A)
Aran=cm_net(sum(A>0));
    
    % Reweight edges
    Dseq=sum(A)';
    for iter=1:50,
        Aran_prev=Aran;
        for iter_node=1:size(A,1),
            tmpsum=sum(Aran(iter_node,:));
            Aran(iter_node,:)=Aran(iter_node,:)/tmpsum*Dseq(iter_node);
            Aran(:,iter_node)=Aran(iter_node,:)';
        end;
        %fprintf('Total change iteration %d: %f\n',iter,sum(sum(abs(Aran-Aran_prev))));
    end;
end