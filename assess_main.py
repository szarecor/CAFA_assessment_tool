#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
#from precrec.precRec import PrecREC, read_benchmark
from precrec.precRec import read_benchmark
from precrec.precRec import PrecREC
#from precision_recall_calculator import PrecisionRecallCalculator as PrecREC
#from precrec.GOPred import GOPred
from go_prediction.go_prediction import GeneOntologyPrediction as GOPred
import os
import sys
import errno
import gc
import yaml


def get_namespace_index(namespace):
    """
    convert namespace into indices
    """
    num = None
    if namespace == "BPO" or namespace == "bpo":
        num = 0
    elif namespace == "MFO" or namespace == "mfo":
        num = 1
    elif namespace == "CCO" or namespace == "cco":
        num = 2
    else:
        raise ValueError("name space not found, check prediction files")
        print(namespace)
    return num


def convert_taxonomy_name(taxonomy_id: str) -> str:
    # convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
    taxonomy_id = str(taxonomy_id)

    taxonomy_table = {
        "10116": "RAT",
        "9606": "HUMAN",
        "3702": "ARATH",
        "7955": "DANRE",
        "44689": "DICDI",
        "7227": "DROME",
        "83333": "ECOLI",
        "10090": "MOUSE",
        "208963": "PSEAE",
        "237561": "CANAX",
        "559292": "YEAST",
        "284812": "SCHPO",
        "8355": "XENLA",
        "224308": "BACSU",
        "99287": "SALTY",
        "243232": "METJA",
        "321314": "SALCH",
        "160488": "PSEPK",
        "223283": "PSESM",
        "85962": "HELPY",
        "243273": "MYCGE",
        "170187": "STRPN",
        "273057": "SULSO",
        "all": "all",
        "prokarya": "prokarya",
        "eukarya": "eukarya",
    }
    return taxonomy_table.get(taxonomy_id)


def convert_benchmark_type(old_type: str) -> str:
    types = {
        'type1': 'NK',
        'type2': 'LK',
        'all': 'All',
    }
    return types.get(old_type)


def extant_file(x):
    if not os.path.isfile(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    else:
        return open(x, "r")

'''
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
'''

def read_config():
    parser = argparse.ArgumentParser(
        description="Precision- Recall assessment for CAFA predictions.",
    )

    parser.add_argument("config_stream", type=extant_file, help="Configuration file")
    # CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    # example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    # If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)
    args = parser.parse_args()
    # Load config file to dictionary
    try:

        config_dict = yaml.safe_load(args.config_stream).get("assess")
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

    obo_path = config_dict.get("obo")
    benchmark_folder = config_dict.get("benchmark")
    results_folder = config_dict.get("results")
    prediction_file = config_dict.get("file")
    return obo_path, benchmark_folder, results_folder, prediction_file


def main():
    from pathlib import Path
    # Pretty terminal output!:
    import colorful as cf
    cf.use_style("monokai")
    color1 = cf.blue
    color2 = cf.green

    # Read Config
    # TODO: Get the metadata in a less dorky way! Use the code from the cafa-format-check repo
    obo_filepath_str, benchmark_folder_str, results_folder_str, prediction_filepath_str = read_config()
    team, model, taxonomy = prediction_filepath_str.rstrip(".txt").split("_")
    keywords = ["sequence alignment"]
    taxonomy_str = convert_taxonomy_name(taxonomy)

    Path(f"{results_folder_str}/precision_recall/").mkdir(parents=True, exist_ok=True)

    predictions_handle = open(prediction_filepath_str, "r")

    bio_process_predictions_path = Path(
        f"{results_folder_str}/{team}_{model}_{taxonomy}_BPO.txt"
    )
    cellular_component_predictions_path = Path(
        f"{results_folder_str}/{team}_{model}_{taxonomy}_CCO.txt"
    )
    molecular_function_predictions_path = Path(
        f"{results_folder_str}/{team}_{model}_{taxonomy}_MFO.txt"
    )

    # Before doing the work of splitting the raw predictions into separate files based on ontologies,
    # check to be sure that it hasn't already been done:
    # TODO: This might not be valid! A team might not have predictions for all three ontologies, correct?
    if not all((bio_process_predictions_path.exists(), cellular_component_predictions_path.exists(), molecular_function_predictions_path.exists())):
        all_predictions = GOPred()
        split_predictions = all_predictions.split_predictions_by_namespace(obo_filepath_str, predictions_handle)
        # TODO: does it make sense to write these separate ontology files to disk?
        with bio_process_predictions_path.open(mode="w") as bio_process_handle:
            bpo_predictions = split_predictions.get("biological process")
            bio_process_handle.write(
                "\n".join([str(prediction) for prediction in bpo_predictions])
            )

        with cellular_component_predictions_path.open(mode="w") as cellular_component_handle:
            cco_predictions = split_predictions.get("cellular component")
            cellular_component_handle.write(
                "\n".join([str(prediction) for prediction in cco_predictions])
            )
        with molecular_function_predictions_path.open(mode="w") as molecular_function_handle:
            mfo_predictions = split_predictions.get("molecular function")
            molecular_function_handle.write(
                "\n".join([str(prediction) for prediction in mfo_predictions])
            )

    print()
    print(" EVALUATING:", color1(prediction_filepath_str))
    print("AUTHOR/TEAM:", color1(team))
    print("      MODEL:", color1(model))
    print("   KEYWORDS:", color1(", ".join(keywords)))
    print("    SPECIES:", color1(f"{taxonomy} ({taxonomy_str})"))

    # Make some new files to hold the evaluation results:
    base_results_filename = prediction_filepath_str.rstrip(".txt")
    results_path_str = f"{results_folder_str}/{base_results_filename}_results.txt"
    precision_recall_path_str = f"{results_folder_str}/precision_recall/{base_results_filename}_precision_recall.txt"

    with open(results_path_str, "w") as results_handle, open(precision_recall_path_str, "w") as precision_recall_handle:

        results_handle.write(f"AUTHOR/TEAM: {team}")
        results_handle.write(f"MODEL: {model}")
        results_handle.write(f"KEYWORDS: {keywords}")
        results_handle.write(f"Species: {taxonomy}")
        results_handle.write("Ontology\tType\tMode\t | Fmax\tThreshold\tCoverage\n")

        indent = " " * 3

        for ontology in ("bpo", "cco", "mfo"):
            ontology_path = f"{results_folder_str}/{team}_{model}_{taxonomy}_{ontology.upper()}.txt"
            print(color2("============================================\n"))
            print("ONTOLOGY:", color1(ontology.upper()))
            for benchmark_type in ("type1", "type2"):

                benchmark_type_str = convert_benchmark_type(benchmark_type)

                print(f"{indent}BENCHMARK TYPE:", color1(f"{benchmark_type_str} ({benchmark_type})"))

                benchmark, obo_count_dict = read_benchmark(
                    namespace=ontology,
                    species=taxonomy_str,
                    types=benchmark_type,
                    fullbenchmarkfolder=benchmark_folder_str,
                    obopath=obo_filepath_str
                )
                if benchmark is None:
                    sys.stderr.write("No benchmark is available for the input species and type")

                c = PrecREC(
                    benchmark=benchmark,
                    #ontology_specific_prediction_path=ontology_path,
                    #ontology_term_count=obo_count_dict.get(ontology)
                    os_pred_path=ontology_path,
                    obocounts=obo_count_dict[ontology]
                )

                if not c.exists:
                    continue

                for mode in ("partial", "full"):
                    print(indent*2, "MODE:", color1(mode))
                    precision, recall, fmax, threshold = c.Fmax_output(mode)
                    coverage = c.coverage()

                    print(indent*3, "      FMAX:", color1(fmax))
                    print(indent*3, " THRESHOLD:", color1(threshold))
                    print(indent*3, "  COVERAGE:", color1(coverage))
                    #print('{:>20}'.format(cf.green(threshold)))
                    print("")

                    results_handle.write(
                        f"{ontology}\t{benchmark_type_str}\t{mode}\t | {fmax}\t{threshold}\t{coverage}\n"
                    )
                    precision_recall_handle.write(
                        f">{ontology}\t{benchmark_type_str}\t{mode}\n |"
                    )
                    precision_recall_handle.write(" ".join([str(i) for i in precision]) + "\n")
                    precision_recall_handle.write(" ".join([str(i) for i in recall]) + "\n")


if __name__ == "__main__":
    main()
