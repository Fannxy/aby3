# clear the records.
rm -r ./aby3-GraphQuery/record_offline/*;

# clear the results.
rm -r ./privGraph/results_offline/*;

python ./privGraph/graph_format_integration_benchmark.py
python ./privGraph/graph_format_record_analysis.py --target privGraph --stage offline

python ./privGraph/graph_format_integration_benchmark.py --target adjmat
python ./privGraph/graph_format_record_analysis.py --target adjmat --stage offline

python ./privGraph/graph_format_integration_benchmark.py --target edgelist
python ./privGraph/graph_format_record_analysis.py --target edgelist --stage offline