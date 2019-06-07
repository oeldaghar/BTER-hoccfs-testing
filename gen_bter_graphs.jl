#generating bter graphs

include("bter_higher_order_testing.jl")
using JLD, HDF5, FileIO

loc = "C:\\Users\\Omar\\Desktop\\Higher-Order_Clustering_Graphs_Testing"
cd(loc)

loc1 = "C:\\Users\\Omar\\Desktop\\Higher-Order_Clustering_Graphs_Testing\\bter_graphs"


ngraphs = 100 #num of bter graphs to generate
for fname in readdir()
  println(fname)
  if fname[end-9:end]==".txt-stats"
    #load data and measure stats for bter
    data = load(fname)["data"]

    maxd = maximum(data[:,1])
    deg_dist = fit(Histogram,data[:,1],nbins=maxd+1)
    ccpd = ccd(data[:,1],data[:,2])
    for j = 1:ngraphs
      println(j)
      E = bter(deg_dist.weights,ccpd)

      #removing self loops
      #TODO double check what they did. maybe increase one of the indices by 1 to avoid self loops
      inds = []
      for i=1:size(E)[1]
        if E[i,1]==E[i,2]
          push!(inds,i)
        end
      end

      E = E[setdiff(1:size(E)[1],inds),:]

    #save edge list
    #TODO this is in the wrong dir. either keep changing directories or use full path name above
      save(File(format"JLD","$(fname[1:end-10])-bter-$(j).txt"),"edges",E)
    end
  end
end
