import inferelator
from inferelator import inferelator_workflow


worker = inferelator_workflow(regression="bbsr", workflow="tfa")
worker.set_file_paths(input_dir="/Users/koenvandenberge/Dropbox/research/GRN/evaluateGRN/yeastGRN/data/",
                    output_dir="/Users/koenvandenberge/Dropbox/research/GRN/evaluateGRN/yeastGRN/scripts/inferNetwork/output_inferelator_bbsr",
                    expression_matrix_file="103118_SS_Data.tsv.gz",
                        tf_names_file="tf_names_yeastract.txt",
                        priors_file="YEASTRACT_20190713_BOTH.tsv",
                        gold_standard_file="YEASTRACT_20190713_BOTH.tsv")
worker.set_run_parameters(num_bootstraps=5, random_seed=42)
network_result = worker.run()

network_result.network.to_csv('/Users/koenvandenberge/Dropbox/research/GRN/evaluateGRN/yeastGRN/scripts/inferNetwork/output_inferelator_bbsr/networkResult.csv')

