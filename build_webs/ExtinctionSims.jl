
using CSV
using DataFrames
using Extinctions
using JLD2
using SpeciesInteractionNetworks

networks = load_object("networks.jlds")

extinctions = load_object("pfim_metaweb.jld2")
filter!(:n_rep => x -> x == 1, extinctions)

# get edgelists
# Filter for G2 (post-extinction) to serve as the "ground truth" for comparison
filter!(:time => x -> x == "G1", networks)# Filter for G2 (post-extinction) to serve as the "ground truth" for comparison
filter!(:model => x -> x == "pfim_metaweb", networks)

#N_post = networks.network[1]
#N_pre = networks.network[1]

#post_edgelist = interactions(networks.network[1])
#pre_edgelist = interactions(networks.network[1])

# Identify species overlap between the two time periods
shared_spp = intersect(species(N_post), species(N_pre))
# The "target richness" used to find the point in the simulation that matches the real world
post_rich = length(shared_spp)

edgelist_df = DataFrame(extinction_mechanism=Any[], consumer=Any[], resource=Any[]);

for i in 1:13

    _ext_seq = extinctions.extinction_seq[i]

    if length(_ext_seq) <= 1
        continue
    end

    net_ind = findfirst(x -> x < post_rich, richness.(_ext_seq)) - 1

    if !(net_ind isa Int) || richness(_ext_seq[net_ind]) == 0
        continue
    end

    sim_N = _ext_seq[net_ind]

    intxns = interactions(sim_N)

    d = DataFrame(
        extinction_mechanism=extinctions.extinction_mechanism[i],
        consumer=getindex.(intxns, 1),
        resource=getindex.(intxns, 2)
    )

    append!(edgelist_df, d)

end

CSV.write("extinction_edgelists.csv",
    edgelist_df)


df_post = DataFrame(
    stage="post",
    consumer=getindex.(post_edgelist, 1),
    resource=getindex.(post_edgelist, 2)
)

df_pre = DataFrame(
    stage="pre",
    consumer=getindex.(pre_edgelist, 1),
    resource=getindex.(pre_edgelist, 2)
)

CSV.write("post_edgelists.csv",
    df_post)

CSV.write("pre_edgelists.csv",
    df_pre)