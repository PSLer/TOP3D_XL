%% Extract diagonal entries (D) of stiffness matrix (K) at Level-l
D_l=zeros(numNodes_l,3); 
eNodMatTmp=eNodMat_l'; eNodMatTmp=eNodMatTmp(:);
if 1==l
diagKe = diag(Ke0);
diagKe = diagKe(:).*E(:)'; %%E: Young's Modulus
eleWithFixedDOFs=find(mapUniqueKes_>0);
eleWithFixedDOFsLocal=mapUniqueKes_(eleWithFixedDOFs);
for kk=1:numel(eleWithFixedDOFs)
kKeFreeDOFs=reshape(uniqueKesFree_(:,eleWithFixedDOFsLocal(kk)),24,24);
kKeFixedDOFs=reshape(uniqueKesFixed_(:,eleWithFixedDOFsLocal(kk)),24,24);
diagKe(:,eleWithFixedDOFs(kk))=diag(kKeFreeDOFs)*E(eleWithFixedDOFs(kk))+diag(kKeFixedDOFs);
end
tmp=diagKe(1:3:end,:); tmp=tmp(:);
D_l(:,1)=D_l(:,1)+accumarray(eNodMatTmp,tmp,[numNodes_l,1]);
tmp=diagKe(2:3:end,:);tmp=tmp(:);
D_l(:,2)=D_l(:,2)+accumarray(eNodMatTmp,tmp,[numNodes_l,1]);
tmp=diagKe(3:3:end,:);tmp=tmp(:);
D_l(:,3)=D_l(:,3)+accumarray(eNodMatTmp,tmp,[numNodes_l,1]);
else
%% 'Ks' stores element stiffness matrices on Level-l (L>1) in the layout 24 x 24 x numElements_l
KsTmp=reshape(Ks,24*24,numElements_l);
diagKeBlock=KsTmp(1:25:(24*24),:);
tmp=diagKeBlock(1:3:end,:); tmp=tmp(:);
D_l(:,1)=D_l(:,1)+accumarray(eNodMatTmp,tmp,[numNodes_l,1]);
tmp = diagKeBlock(2:3:end,:); tmp=tmp(:);
D_l(:,2)=D_l(:,2)+accumarray(eNodMatTmp,tmp,[numNodes_l,1]);
tmp=diagKeBlock(3:3:end,:); tmp=tmp(:);
D_l(:,3)=D_l(:,3)+accumarray(eNodMatTmp,tmp,[numNodes_l,1]);
end
D_l=reshape(D_l',numDOFs_l,1);