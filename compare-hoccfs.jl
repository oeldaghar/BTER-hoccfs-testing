using JLD,HDF5,FileIO
using StatsBase,Plots
pyplot()

#compare hoccf distributions
function hoccpd(degs,ccfs)
    maxd = Int(maximum(degs))
    ccpd = zeros(Float64,maxd)

    for i=1:maxd
        inds = findall(degs.==i)
        if length(inds)>0
            ccpd[i] = sum(ccfs[inds])/length(inds)
        end
    end
    return ccpd
end

function compare_deg_dist(degs,degs1)
    maxd = maximum(degs)
    maxd1 = maximum(degs1)

    dd = fit(Histogram,degs,nbins=maxd+1)
    dd1 = fit(Histogram,degs1,nbins=maxd1+1)

    #remove empty degs
    inds = findall((dd.weights).==0)
    inds1 = findall((dd1.weights).==0)

    edges = collect(dd.edges[1])[1:end-1][Int.(setdiff(1:maxd,inds))]
    weights = collect(dd.weights)[Int.(setdiff(1:maxd,inds))]

    edges1 = collect(dd1.edges[1])[1:end-1][Int.(setdiff(1:maxd1,inds1))]
    weights1 = collect(dd1.weights)[Int.(setdiff(1:maxd1,inds1))]

    Plots.scatter(edges,weights)
    Plots.scatter!(edges1,weights1,c = :red)
    gui()
end

function compare_hoccpd(degs,degs1,ccpd,ccpd1)
    maxd = maximum(degs)
    maxd1 = maximum(degs1)

    dd = fit(Histogram,degs,nbins=maxd+1)
    dd1 = fit(Histogram,degs1,nbins=maxd1+1)

    #remove empty degs
    inds = findall((dd.weights).==0)
    inds1 = findall((dd1.weights).==0)

    ccpd = ccpd[setdiff(1:maxd,inds)]
    ccpd1 = ccpd1[setdiff(1:maxd1,inds1)]

    Plots.scatter(setdiff(1:maxd,inds),ccpd)
    Plots.scatter!(setdiff(1:maxd1,inds1),ccpd1,c=:red)
end

function hoccfspd(data)
    degs = data[:,1]
    ccfs3 = data[:,2]
    ccfs4 = data[:,3]
    ccfs5 = data[:,4]
    maxd = Int(maximum(degs))

    result = zeros(maxd,4)
    for i=1:maxd
        inds = findall(degs.==i)
        n = length(inds)
        if n>0
            result[i,1] = n
            result[i,2] = sum(ccfs3[inds])/n
            result[i,3] = sum(ccfs4[inds])/n
            result[i,4] = sum(ccfs5[inds])/n
        end
    end
    return result
end

function remove_empty_degs(hoccfspd_data)
    maxd = size(hoccfspd_data)[1]
    inds = []
    for i=1:maxd
        if hoccfspd_data[i,1]==0
            push!(inds,i)
        end
    end
    return setdiff(1:maxd,inds),hoccfspd_data[setdiff(1:maxd,inds),:]
end

#loading real data
loc = "C:\\Users\\Omar\\Desktop\\Higher-Order_Clustering_Graphs_Testing\\processed_graphs"
bter_loc = "C:\\Users\\Omar\\Desktop\\Higher-Order_Clustering_Graphs_Testing\\bter_graphs"
their_bter_loc = "C:\\Users\\Omar\\Desktop\\Higher-Order_Clustering_Graphs_Testing\\their-bter-graphs"
cd(loc)
gnames = readdir()[findall(endswith.(readdir(),".txt-stats"))]
#fname = gnames[3]
for i = 1:length(gnames)
    cd(loc)
    fname = gnames[i]#graph
    data = load(fname)["data"] #real data
    cd(bter_loc)
    bter_files = filter(x->startswith(x,fname[1:end-10]) && endswith(x,".txt-stats"),readdir())
    bter_data = load(bter_files[rand(1:length(bter_files))])["data"]#rand bter graph from directory
    cd(their_bter_loc)
    their_bter_files = filter(x->startswith(x,"bter-$(fname[1:end-10])") && endswith(x,".txt-stats"),readdir())
    their_bter_data = load(their_bter_files[rand(1:length(their_bter_files))])["data"]
    x1,y1 = remove_empty_degs(hoccfspd(data))
    x2,y2 = remove_empty_degs(hoccfspd(bter_data))
    x3,y3 = remove_empty_degs(hoccfspd(their_bter_data))
    for k=2:4
        scatter(x1,y1[:,k],label="original")
        scatter!(x2,y2[:,k],c=:red,label="my bter sim")
        scatter!(x3,y3[:,k],c=:green,label="their bter sim")
        xlabel!("degree")
        ylabel!("avg HOCCF per deg")
        title!("HOCCF per degree dist for k=$(k)")
        gui()
    end
end

fname = gnames[1]
cd(loc)
data = load(fname)["data"] #real data
cd(bter_loc)
bter_files = filter(x->startswith(x,fname[1:end-10]) && endswith(x,".txt-stats"),readdir())
bter_data = load(bter_files[rand(1:length(bter_files))])["data"]#rand bter graph from directory
cd(their_bter_loc)
their_bter_files = filter(x->startswith(x,"bter-$(fname[1:end-10])") && endswith(x,".txt-stats"),readdir())
their_bter_data = load(their_bter_files[rand(1:length(their_bter_files))])["data"]
x1,y1 = remove_empty_degs(hoccfspd(data))
x2,y2 = remove_empty_degs(hoccfspd(bter_data))
x3,y3 = remove_empty_degs(hoccfspd(their_bter_data))
k = 4
scatter(x1,y1[:,k],label="original")
scatter!(x2,y2[:,k],c=:red,label="my bter sim")
scatter!(x3,y3[:,k],c=:black,label="their bter sim")
xlabel!("degree")
ylabel!("avg HOCCF per deg")
title!("$(fname[1:end-10]) HOCCF per degree dist for k=$(k)")
gui()



cd(loc)
(A,inds) = read_undir_graph_txt(fname[1:end-6],true)
#data = load(fname)["data"]
degs = A*ones(size(A)[1])
ccfs = clustercoeffs(A)
ccd(degs,ccfs)
