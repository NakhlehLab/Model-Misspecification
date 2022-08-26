using PhyloNetworks
using QuartetNetworkGoodnessFit, DataFrames, CSV
using DelimitedFiles

function computePvalue1replicate(file_path::AbstractString, candidate_net_string::AbstractString, optimize::Bool, output_df_path::AbstractString, output_res_path::AbstractString, ingroup::Bool)
	qCF = DataFrame(CSV.File(file_path), copycols=false);
	net0 = readTopology(candidate_net_string);
	if ingroup
		deleteleaf!(net0, "Z");
	end

	if optimize
		res0 = quarnetGoFtest!(net0, qCF, true; seed=201, nsim=5);
		net0 = res0[5]
	end
	res0 = quarnetGoFtest!(net0, qCF, false; seed=234, nsim=200);
	CSV.write(output_df_path, qCF)
	open(output_res_path, "w") do file
		writedlm(file, res0[[1,2,3]], ",")
	end
	return res0
end

function read_ML_nettopo(ML_path::AbstractString)
	open(ML_path) do f
        line = 1
        find=false
        res=""
        while !eof(f)
            x = readline(f)
            if occursin("Inferred Network #1", x)
                find=true
            elseif find
                res=x
                break
            end
            line += 1
        end
        return res
    end
end

function computeExpCF(df_obsCF::DataFrame, output_path::AbstractString)
	dcf = readTableCF(df_obsCF)
	res = quarnetGoFtest!(net0, dcf, true; seed=201, nsim=3);
	df_out=writeExpCF(dcf.quartet)
	CSV.write(output_path, df_out)
end


function main1_args()
	print(ARGS)
	path=ARGS[1]
	net_path=ARGS[2]
	output_df_path=ARGS[3]
	output_res_path=ARGS[4]
	ingroup=(lowercase(ARGS[5]) == "true")
	net_string = read_ML_nettopo(net_path)
	res=computePvalue1replicate(path, net_string, true, output_df_path, output_res_path, ingroup)
	print(res)
end


main1_args()