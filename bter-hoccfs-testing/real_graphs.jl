using HigherOrderClustering
using JLD, HDF5, FileIO

loc = "C:\\Users\\Omar\\Desktop\\Higher-Order_Clustering_Graphs_Testing\\processed_graphs"
cd(loc)

#measuring data from true graphs
for file in readdir()
   if file[end-3:end]==".txt"
      println(file)
      (A,inds) = read_undir_graph_txt(file, true)
      nnodes = size(A)[1]
      degs = A*ones(nnodes)
      ccfs2 = higher_order_ccfs(A,2) #classical local clustering 
      ccfs3 = higher_order_ccfs(A,3) 
      ccfs4 = higher_order_ccfs(A,4)  	   

      #storing and saving data
      data = zeros(Float64,nnodes,4)
      data[:,1] = degs
      data[:,2] = ccfs2.local_hoccfs
      data[:,3] = ccfs3.local_hoccfs
      data[:,4] = ccfs4.local_hoccfs

      save(File(format"JLD","$(file)-stats"),"data",data)
   end
end 
