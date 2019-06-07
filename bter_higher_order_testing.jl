#TODO cite bter stuff and reproduce(?) their copyright and stuff

using DelimitedFiles #for loading data
using StatsBase
using SparseArrays

function bter_setup(nd,cd,beta=1) #nd deg dist #cd ccpd
    #@assert() TODO  check that this is a valid deg seq
    dmax=length(nd)
    #initialize stuff
    id=zeros(dmax)
    wd=zeros(dmax)
    rdfill=zeros(dmax)
    ndfill=zeros(dmax)
    wg=zeros(dmax)
    ig=zeros(dmax)
    bg=zeros(dmax)
    ng=zeros(dmax)
    ndprime=zeros(dmax) #nd'

    #relabel nodes in increasing degree w/ deg 1 nodes at end
    tmp=cumsum(nd[2:end])
    id[2]=1
    id[3:end] = tmp[1:end-1].+1
    id[1] = tmp[end]+1

    #computing nodes of degree greater than d
    tmp = reverse(cumsum(reverse(nd)))
    ndprime[2:end-1] = tmp[3:end]

    #handle degree 1 nodes
    ndfill[1]=nd[1]*beta
    wd[1]=0.5*nd[1]
    rdfill[1]=1

    #main loop
    g = 0
    nfillblk = 0
    intdeg = 0
    for d = 2:dmax

        if nfillblk > 0
            ndfill[d] = min( nfillblk, nd[d] )
            nfillblk = nfillblk - ndfill[d]
            wdfilltmp = 0.5 * ndfill[d] * (d - intdeg)
        else
            ndfill[d] = 0
            wdfilltmp = 0
        end

        ndbulktmp = nd[d] - ndfill[d]

        if ndbulktmp > 0
            g = g + 1
            ig[g] = id[d] + ndfill[d]
            bg[g] = ceil(ndbulktmp / (d+1))
            ng[g] = d+1
            if (bg[g]*(d+1)) > (ndprime[d] + ndbulktmp)
                if bg[g] != 1
                    error("Last group has more than 1 block")
                end
                ng[g] = ndprime[d] + ndbulktmp
            end
            rho = (cd[d])^(1/3)
            intdeg = (ng[g] - 1) * rho
            wdbulktmp = 0.5 * ndbulktmp * (d - intdeg);
            wg[g] = bg[g] * 0.5 * ng[g] * (ng[g] - 1) * log(1/(1-rho));
            nfillblk = bg[g] * ng[g] - ndbulktmp
        else
            wdbulktmp = 0
        end

        wd[d] = wdbulktmp + wdfilltmp
        if wd[d] > 0
            rdfill[d] = wdfilltmp / wd[d]
        else
            rdfill[d] = 0
        end
    end

    #Shorten the group arrays
    ig = ig[1:g]
    wg = wg[1:g]
    bg = bg[1:g]
    ng = ng[1:g]

    return [id,wd,ndfill,rdfill,ig,wg,bg,ng]
end

function random_sample(cnts,nsamples)
    #if !@isdefined(nsamples) #TODO get nsamples parameter working apart from default
    #    nsamples = round(sum(cnts))
    #else
    #    cnts = cnts.*(nsamples/sum(cnts))
    #end
    cnts = cnts.*(nsamples/sum(cnts))
    cumdist = [0; cumsum(cnts)]
    bins = cumdist/cumdist[end]
    #testval = abs(bins[end] - 1)
    #if  testval > eps()
    #    @warn("Last entry of bins is not exactly 1. Diff = %e.", testval)
    #end

    #naive way of doing this. need to check that this is the same
#    h = fit(Histogram,rand(nsamples),bins)
#    result=[]
#    for ind in length(h.weights)
#        if h.weights[ind]>0
#            append!(result,ones(Int64,h.weights[ind])*cnts[ind])
#        end
#    end
#    return result
    #or manually
    result=zeros(Int64,nsamples)

    for ind=1:nsamples
        r=rand(1)[1]
        index=findfirst(bins.>= r)
        result[ind]=index-1
    end
    return result
end


function bter(nd,cd,beta=1)
    id,wd,ndfill,rdfill,ig,wg,bg,ng = bter_setup(nd, cd, beta)

    w1 = sum(wg)
    w2 = sum(wd)
    w = w1+w2
    nsmp = round(w)
    r = rand(Int(nsmp),1)
    s1 = sum(r .< w1/w)
    s2 = nsmp - s1

    #phase 1 samples
    grp_smp = random_sample(wg, s1)
    blk_r = rand(s1,1)
    blk_b = bg[grp_smp]
    blk_i = ig[grp_smp]
    blk_n = ng[grp_smp]
    e1_shift = blk_i + floor.(blk_r .* blk_b) .* blk_n
    e1_r = rand(s1,2)
    e1 = zeros(s1,2)
    e1[:,1] = floor.(e1_r[:,1] .* blk_n) .+ e1_shift
    e1[:,2] = floor.(e1_r[:,2] .* (blk_n .- 1)) .+ e1_shift
    e1[:,2] = e1[:,2] .+ Float64.((e1[:,2] >= e1[:,1]))


    #phase 2 samples
    # Setup
    idfill = id;
    idbulk = id + ndfill;
    ndbulk = nd - ndfill;
    ndbulk[1] = 0;
    # Sample
    d_smp = random_sample(wd, Int(2*s2));
    d_smp = reshape(d_smp, Int(s2), 2);
    r = rand(Int(s2),2);
    tf_fill = r .< rdfill[d_smp];
    e2_shift_fill = idfill[d_smp];
    e2_sz_fill = ndfill[d_smp];
    e2_shift_bulk = idbulk[d_smp];
    e2_sz_bulk = ndbulk[d_smp];

    r = rand(Int(s2),2);
    e2_fill = e2_shift_fill + floor.(r .* e2_sz_fill);
    e2_bulk = e2_shift_bulk + floor.(r .* e2_sz_bulk);
    e2 = tf_fill .* (e2_fill) + (1 .- tf_fill) .* (e2_bulk);

    #TODO combine e1 and e2 into final edge list and make CSC mat
    return [e1; e2]
end

#this function comes from Austin Benson's code for
function read_undir_graph_txt(filename::AbstractString, oneindex::Bool=false)
    # index mapping
    index_map = Dict{Int64,Int64}()
    index_map_vec = Int64[]
    function get_mapped_index(x::Int64)
        if !haskey(index_map, x)
            next_index = length(index_map) + 1
            index_map[x] = next_index
            push!(index_map_vec, x)
            return next_index
        end
        return index_map[x]
    end

    # Read data
    I = Int64[]
    J = Int64[]
    open(filename) do f
        for line in eachline(f)
            # Skip lines starting with '#' or '%'
            if line[1] == '#' || line[1] == '%'; continue; end
            edge = split(line)
            u = parse(Int64, edge[1])
            v = parse(Int64, edge[2])
            if oneindex
                push!(I, get_mapped_index(u))
                push!(J, get_mapped_index(v))
            else
                push!(I, u)
                push!(J, v)
            end
        end
    end

    # Form adjacency matrix
    if min(minimum(I), minimum(J)) < 1
        error("Minimum node value is less than 1. Try setting ondeindex parameter to true.")
    end
    n = max(maximum(I), maximum(J))
    A = convert(SparseMatrixCSC{Int64,Int64}, sparse(I, J, ones(length(I)), n, n))
    A = max.(A, A')

    if oneindex; return (A, index_map_vec); end
    return (A, collect(1:n))
end


function ccd(degs,local_ccfs)
    #degs: degs for each node, local_ccfs: local_ccfs for each node
    maxd = Int(maximum(degs))
    ccd = zeros(maxd)

    for i = 1:maxd
        inds = (degs .== i) #nodes of deg i
        ccd[i] = sum(local_ccfs[inds])/length(inds) # avg local cc for nodes of deg i
    end
    return ccd
end

#=
#TODO there is a bug with Austin Benson's code for reading in a graph. see CA-GrQc.
degree sum should be 28992 but is 28980 even after processing
=#
