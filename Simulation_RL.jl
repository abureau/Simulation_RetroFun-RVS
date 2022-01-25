#Julia Script to parse and convert vcf file to genotypes
using(Pkg)
using(GeneticVariation)
using(DataFrames)
using(DelimitedFiles)
using(Distributions)
using(RData)
using(CSV)
using(Combinatorics)
using(LinearAlgebra)
using(StatsBase)
using(Plots)
using(Glob)
using(Random)

Random.seed!(1234)

function read_vcf(path, sep='\t', comment_char='#')
    #This function loads a vcf file and format it for further manipulation
    vcf_header = readlines(path)[2]
    columns = split(strip(vcf_header,[comment_char]),sep)
    vcf_file = DataFrame(readdlm(path, sep, Any, '\n', comments=true, comment_char=comment_char), :auto)
    rename!(vcf_file,columns)
end

function convert_geno(df_vcf)
    #This function converts "1|0", "0|1","0|0" or "1|1" into genotypes for a given risk genetic model

    ind = size(df_vcf,2)
    df_vcf[:,10:ind].=ifelse.(df_vcf[:,10:ind].==("1|0"), "0|1", df_vcf[:,10:ind])


    df_vcf[:,10:ind].=ifelse.(df_vcf[:,10:ind].==("0|1"), 1, df_vcf[:,10:ind])
    df_vcf[:,10:ind].=ifelse.(df_vcf[:,10:ind].==("1|1"), 2, df_vcf[:,10:ind])
    df_vcf[:,10:ind].=ifelse.(df_vcf[:,10:ind].==("0|0"), 0, df_vcf[:,10:ind])


end

function read_ped(path, sep=' ')
    ped_file = DataFrame(readdlm(path, sep, Any, '\n', comments=true, comment_char='#'), :auto)

    n_col = size(ped_file, 2) - 6

    ped_info_columns = Vector{String}(["FamID", "ID", "DadID", "MomID", "sex", "status"])
    locus_name = map(string,repeat(["Locus"], n_col),[1:1:n_col;])

    rename!(ped_file, vcat(ped_info_columns, locus_name))
end

function convert_geno_from_ped(pedfile)

    locus_col = size(pedfile,2)
    pedfile[:,7:locus_col].=ifelse.(pedfile[:,7:locus_col].==(1), 0, pedfile[:,7:locus_col])
    pedfile[:,7:locus_col].=ifelse.(pedfile[:,7:locus_col].==(2), 1, pedfile[:,7:locus_col])
    pedfile[:,7:locus_col].=ifelse.(pedfile[:,7:locus_col].==(3), 2, pedfile[:,7:locus_col])

    hcat(pedfile[:, 1:6],  pedfile[:,7:locus_col])
end

function combine_genos(pedfile)
    ncol = size(pedfile,2)
    z = zip([7:2:ncol;],[8:2:ncol;])
    ped_file_agg = DataFrame()

    for (k,(i,j)) in enumerate(z)
        insertcols!(ped_file_agg,k,:"Locus".*string(k)=>ped_file[:,i].+ped_file[:,j])
    end

    ped_file_agg = hcat(ped_file[:,Between("FamID","status")], ped_file_agg)

end

function aggregate_by_fam(pedfile)

    pedfile_affected = filter(:status => ==(2), pedfile)

    ncol = size(pedfile, 2)
    ncol_genotypes = ncol - 6

    locus_cols = map(string,repeat(["Locus"], ncol_genotypes),[1:1:ncol_genotypes;])

    pedfile_affected_agg = combine(groupby(pedfile_affected[:, [1:2;7:ncol]], ["FamID"]),  locus_cols .=> sum)
    agg_by_fam = hcat(pedfile_affected_agg[:,"FamID"],sum.(eachrow(pedfile_affected_agg[:,Not("FamID")])))

    Dict("Agg"=> agg_by_fam, "Row"=> pedfile_affected_agg)
end

function r_zero_to_NaN(x)
    replace(x, 0=>NaN)
end

function r_NaN_to_zero(x)
    replace(x, NaN=>0)
end

function Variance_Stats_burden(pedfile)
    sharing_prob = Dict()
    var_genotypes = Dict()
    cov_genotypes = Dict()

    agg_by_fam = aggregate_by_fam(pedfile)

    all_agg_by_fam = get(agg_by_fam,"Agg", 0)
    row_genos_by_fam = get(agg_by_fam,"Row", 0)

    row_genos_by_fam.FamID = string.(row_genos_by_fam.FamID)
    row_genos_by_fam.FamID = ifelse.(isuppercase.(first.(row_genos_by_fam.FamID)), chop.(row_genos_by_fam.FamID, tail=0, head=1), row_genos_by_fam.FamID)

    non_null_fam = string.(all_agg_by_fam[all_agg_by_fam[:,2] .>0,:1])
    fam_in_sim = Set([ifelse(isuppercase(first(famid)) == (true),chop(famid,head=1,tail=0), famid) for famid in non_null_fam ])

    #agg_by_fam_non_null_genos = filter(:FamID => in(non_null_fam),row_genos_by_fam)


    null_value_for_fam_in_Sim = Dict(sub_key => null_value[sub_key] for sub_key in fam_in_sim)

    row_genos_by_fam
end

cd("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\")

SPAP_prob = load("SPAPsimpleprob.RData", convert=true)
prob_sharing = get(SPAP_prob,"forsim.pattern.prob.list",0)
N_pattern = get(SPAP_prob,"forsim.N.list",0)

MAF_var = CSV.File("weights_subset_haplo_sfs_2_rare.txt", header=0) |> DataFrame

weights_var = diagm(pdf.(Beta(1,25),diag(diagm(MAF_var[:,:1]))))

null_value = Dict(key.=>sum(prob_sharing[key].*N_pattern[key]) for key in keys(prob_sharing))

var_alleles = Dict()
cov_alleles = Dict()

for key in keys(prob_sharing)
    push!(var_alleles,key => sum(prob_sharing[key].*N_pattern[key].^2)-null_value[key].^2)

    joint = [minimum([i,j]) for j in prob_sharing[key], i in prob_sharing[key]]/sum([minimum([i,j]) for j in prob_sharing[key], i in prob_sharing[key]])
    marginal = sum.(eachcol(joint))
    moy = sum(N_pattern[key].*marginal)

    joint_config = [[i*j] for j in N_pattern[key], i in N_pattern[key]]
    joint_config_joint = sum(joint_config.*joint) .- moy^2

    push!(cov_alleles, key => joint_config_joint)
end

cov_alleles_approx = Dict(key => minimum(vcat(cov_alleles[key],var_genotypes[key])) for key in keys(prob_sharing))

Stats = []
Var = []

Array_df_genotypes = []

for f in glob("Null/CRH/rep*.ped")

    ped_file = read_ped("$f", ' ')
    ped_file = convert_geno_from_ped(ped_file)
    ped_file_agg = combine_genos(ped_file)

    ped_file_affected = filter(:status => ==(2), ped_file_agg)
    #push!(Array_df_genotypes, ped_file_affected)

    #locus_cols = map(string,repeat(["Locus"], (size(ped_file)-6)/2),[1:1:(size(ped_file)-6)/2;])
    #ped_file_affected_agg = combine(groupby(ped_file_affected[:, [1:2;7:size(ped_file)]], ["FamID"]),  locus_cols .=> sum)
    agg_by_fam = aggregate_by_fam(ped_file_affected)

    all_agg_by_fam = get(agg_by_fam,"Agg", 0)
    all_agg_by_fam[:,:1] = string.(all_agg_by_fam[:,:1])
    all_agg_by_fam[:,:1] = ifelse.(isuppercase.(first.(all_agg_by_fam[:,:1])), chop.(all_agg_by_fam[:,:1], tail=0, head=1), all_agg_by_fam[:,:1])

    #fam_non_null_genos = Set(all_agg_by_fam[all_agg_by_fam[:,:2] .>0,:1])

    row_genos_by_fam = get(agg_by_fam,"Row", 0)
    row_genos_by_fam.FamID = string.(row_genos_by_fam.FamID)
    row_genos_by_fam.FamID = ifelse.(isuppercase.(first.(row_genos_by_fam.FamID)), chop.(row_genos_by_fam.FamID, tail=0, head=1), row_genos_by_fam.FamID)

    #agg_by_fam = hcat(ped_file_affected_agg[:,"FamID"],sum.(eachrow(ped_file_affected_agg[:,Not("FamID")])))
    #fam_non_null_genos = all_agg_by_fam[all_agg_by_fam[:,2] .>0,:1]
    agg_by_fam_non_null_genos = row_genos_by_fam[findall(all_agg_by_fam[:,:2] .>0),:]

    #fam_in_sim = Set([ifelse(isuppercase(first(famid)) == (true),chop(famid,head=1,tail=0), famid) for famid in fam_non_null_genos ])
    #null_value = Dict(key.=>sum(prob_sharing[key].*N_pattern[key]) for key in keys(prob_sharing))
    #null_value_for_fam_in_Sim = Dict(sub_key => null_value[sub_key] for sub_key in fam_in_sim)

    t = agg_by_fam_non_null_genos

    mapcols!(r_zero_to_NaN,t)

    array_dfs_fam = []

    for r in eachrow(t)
        #push!(array_dfs_fam,filter(:FamID => ==(r.FamID), t)[:,2:size(t,2)] .-= get(null_value, r.FamID, 0))
        push!(array_dfs_fam,collect(r[2:size(t,2)]) .- get(null_value, r.FamID, 0))
    end

    #concat_dfs_fam = reduce(vcat, array_dfs_fam)
    concat_dfs_fam = DataFrame(convert(Array,hcat(array_dfs_fam...)'), :auto)

    mapcols!(r_NaN_to_zero, concat_dfs_fam)

    rename!(concat_dfs_fam,[string("SNP", x) for x in 1:1:size(concat_dfs_fam,2)])

    push!(Array_df_genotypes, hcat(t.FamID,concat_dfs_fam))
    push!(Stats, sum(sum.(eachrow(Matrix(concat_dfs_fam)*weights_var))))

    var_alleles_tmp = Dict(fam => [] for fam in keys(null_value))
    covar_alleles_tmp = Dict(fam => [] for fam in keys(null_value))

    for r in eachrow(agg_by_fam_non_null_genos)

        genos = collect(r[2:size(r,1)])[collect(r[2:size(r,1)]) .> 0]
        index_genos_non_null = findall(collect(r[2:size(r,1)]) .> 0)

        subset_weight = diag(weights_var[index_genos_non_null,index_genos_non_null])

        push!(var_alleles_tmp, r.FamID => [sum(subset_weight.^2*get(var_alleles, r.FamID, 0))])


        n_genos = length(genos)
        n_unique_genos = length(unique(genos))

        if n_genos == 1
            push!(covar_alleles_tmp, r.FamID=>[get(covar_alleles_tmp, r.FamID, 0);get(var_alleles_tmp, r.FamID,0)])

        else

            array_sum_weight = []

            for i in combinations(subset_weight,2)
                push!(array_sum_weight,prod(i))
            end

            if n_unique_genos == 1
                push!(covar_alleles_tmp, r.FamID => [get(covar_alleles_tmp, r.FamID, 0);get(var_alleles_tmp, r.FamID, 0).+2*sum(array_sum_weight)*get(var_alleles, r.FamID, 0)])

            else
                push!(covar_alleles_tmp, r.FamID => [get(covar_alleles_tmp, r.FamID, 0);get(var_alleles_tmp, r.FamID, 0).+2*sum(array_sum_weight)*get(cov_alleles_approx, r.FamID,0)])

            end
        end

    end

    push!(Var, sum(vcat(values(covar_alleles_tmp)...)))

end

Stats_with_Var = hcat(Stats,Var)

Stats_from_Parametric_boot = []

for b in 1:1:1000
    re_Q = rand(Normal(0, sqrt(Stats_with_Var[b,:2])), 10000)

    re_Mean = mean(re_Q)
    re_Variance = var(re_Q)
    re_Kurtosis = kurtosis(re_Q)

    re_df=(re_Kurtosis>0)*12/re_Kurtosis+(re_Kurtosis<=0)*100000

    adjusted_Q = (Stats_with_Var[b,:1]-re_Mean)*sqrt(2*re_df)/sqrt(re_Variance)+re_df
    push!(Stats_from_Parametric_boot,cdf(Chisq(re_df),adjusted_Q))

end

mean(Stats_with_Var[:,:1].^2 ./ Stats_with_Var[:,:2])

histogram(Stats_with_Var[:,:1] ./ sqrt.(Stats_with_Var[:,:2]))
histogram(Stats_with_Var[:,:1].^2 ./ Stats_with_Var[:,:2])

pvalues_initial = ccdf(Chisq(1),Stats.^2 ./ Var)

QQPlot = plot(-log10.(collect(1:1:length(pvalues_initial))/length(pvalues_initial)),
     -log10.(sort(pvalues_initial)), seriestype = :scatter,
     )
Plots.abline!(1, 0, line=:dash)
savefig(QQPlot, "C:/Users/loicm/Documents/recherche/Vraisemblance_retrospective/Simulation/output/QQPlot.png")
QQPlot
Plots.abline!(1, 0, line=:dash)
small_pvalues = findall(pvalues_initial.<=0.01)
Stat_nonparam_Bootstrap = Dict(index => [] for index in small_pvalues)

prob_sharing_by_fam = Dict(key => [] for key in keys(prob_sharing))

for famid in keys(prob_sharing)
    for i in 1:1:length(groupinds(groupslices(N_pattern[famid])))
        index = groupinds(groupslices(values(N_pattern[famid])))[i]
        push!(prob_sharing_by_fam, famid => [get(prob_sharing_by_fam, famid, 0);sum(prob_sharing[famid][index])])
    end
    reverse!(prob_sharing_by_fam[famid])
end


for index in small_pvalues

    df = Array_df_genotypes[index]

    fam = df[:,:1]

    n_variants_by_famid = count.(!=(0), eachrow(df[:,2:36]))

    index_variants_by_famid = [parse.(Int8,hcat(split.(String.(x), "SNP")...)[:2,:]) for x in findall.(!=(0), eachrow(df[:,2:36]))]

    for boot in 1:1:10000
        Stat_by_fam = []
        for i in 1:1:length(n_variants_by_famid)
            sample_geno = sample(1:1:length(prob_sharing_by_fam[fam[i]]),Weights(Float64.(prob_sharing_by_fam[fam[i]])), n_variants_by_famid[i], replace=true)

            push!(Stat_by_fam, sum(diag(weights_var[index_variants_by_famid[i], index_variants_by_famid[i]]).*(sample_geno .- null_value[fam[i]])))

        end
        push!(Stat_nonparam_Bootstrap,index => [get(Stat_nonparam_Bootstrap,index,0);sum(Stat_by_fam)])

    end

end


bootstrap_pvalues = [sum(Stat_nonparam_Bootstrap[k].^2 .>= Stats[k]^2)/1000 for k in keys(Stat_nonparam_Bootstrap)]

deleteat!(pvalues_initial, small_pvalues)

pvalues_boot = pvalues_initial[Not(small_pvalues)]

toto = [pvalues_boot;bootstrap_pvalues]

plot(-log10.(collect(1:1:length(toto))/length(toto)),
     -log10.(sort(toto)), seriestype = :scatter,
     )
Plots.abline!(1, 0, line=:dash)

boot_agg_by_fam = aggregate_by_fam(df_boot)

boot_all_agg_by_fam = get(boot_agg_by_fam,"Agg", 0)
boot_all_agg_by_fam[:,:1] = string.(boot_all_agg_by_fam[:,:1])
boot_all_agg_by_fam[:,:1] = ifelse.(isuppercase.(first.(boot_all_agg_by_fam[:,:1])), chop.(boot_all_agg_by_fam[:,:1], tail=0, head=1), boot_all_agg_by_fam[:,:1])

boot_row_genos_by_fam = get(boot_agg_by_fam,"Row", 0)
boot_row_genos_by_fam.FamID = string.(boot_row_genos_by_fam.FamID)
boot_row_genos_by_fam.FamID = ifelse.(isuppercase.(first.(boot_row_genos_by_fam.FamID)), chop.(boot_row_genos_by_fam.FamID, tail=0, head=1), boot_row_genos_by_fam.FamID)

boot_agg_by_fam_non_null_genos = boot_row_genos_by_fam[findall(boot_all_agg_by_fam[:,:2] .>0),:]

boot_t = boot_agg_by_fam_non_null_genos

mapcols!(r_zero_to_NaN,boot_t)

boot_array_dfs_fam = []

for r in eachrow(boot_t)
    #push!(array_dfs_fam,filter(:FamID => ==(r.FamID), t)[:,2:size(t,2)] .-= get(null_value, r.FamID, 0))
    push!(boot_array_dfs_fam,collect(r[2:size(boot_t,2)]) .- get(null_value, r.FamID, 0))
end

#concat_dfs_fam = reduce(vcat, array_dfs_fam)
boot_concat_dfs_fam = DataFrame(convert(Array,hcat(boot_array_dfs_fam...)'), :auto)

mapcols!(r_NaN_to_zero, boot_concat_dfs_fam)
