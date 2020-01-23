#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from precrec.precRec import PrecREC, read_benchmark
#from precrec.GOPred import GOPred
from go_prediction.go_prediction import GOPrediction as GOPred
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


def taxon_name_converter(taxonID):
    # convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
    taxonTable = {
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
    return taxonTable[taxonID]


def typeConverter(oldType):
    if oldType == "type1":
        newType = "NK"
    elif oldType == "type2":
        newType = "LK"
    elif oldType == "all":
        newType = "All"
    return newType


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

        config_dict = yaml.safe_load(args.config_stream)["assess"]
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()

    obo_path = config_dict["obo"]
    benchmark_folder = config_dict["benchmark"]
    results_folder = config_dict["results"]
    prediction_file = config_dict["file"]
    return obo_path, benchmark_folder, results_folder, prediction_file


def main():
    from pathlib import Path
    import colorful as cf

    cf.use_style("monokai")

    # Read Config
    obo_path, benchmark_folder, results_folder, prediction_file = read_config()
    team, model, taxonomy = prediction_file.rstrip(".txt").split("_")
    keywords = ["sequence alignment"]
    taxonomy_str = taxon_name_converter(taxonomy)

    Path(f"{results_folder}/precision_recall/").mkdir(parents=True, exist_ok=True)

    all_predictions = GOPred()
    predictions_handle = open(prediction_file, "r")
    split_predictions = all_predictions.split_predictions_by_namespace(obo_path, predictions_handle)

    # TODO: does it make sense to write these separate ontology files to disk?

    with open(f"{results_folder}/{team}_{model}_{taxonomy}_BPO.txt", "w") as bio_process_handle:
        bpo_predictions = split_predictions.get("biological process")
        bio_process_handle.write(
            # TODO: Do we really need to explicitly called __repr__() ?
            "\n".join([prediction.__repr__() for prediction in bpo_predictions])
        )

    with open(f"{results_folder}/{team}_{model}_{taxonomy}_CCO.txt", "w") as cellular_component_handle:
        cco_predictions = split_predictions.get("cellular component")
        cellular_component_handle.write(
            # TODO: Do we really need to explicitly called __repr__() ?
            "\n".join([prediction.__repr__() for prediction in cco_predictions])
        )

    with open(f"{results_folder}/{team}_{model}_{taxonomy}_MFO.txt", "w") as molecular_function_handle:
        mfo_predictions = split_predictions.get("molecular function")
        molecular_function_handle.write(
            # TODO: Do we really need to explicitly called __repr__() ?
            "\n".join([prediction.__repr__() for prediction in mfo_predictions])
        )

    print()
    print(" EVALUATING:", cf.green(prediction_file))
    print("AUTHOR/TEAM:", cf.green(team))
    print("      MODEL:", cf.green(model))
    print("   KEYWORDS:", cf.green(", ".join(keywords)))
    print("    SPECIES:", cf.green(f"{taxonomy} ({taxonomy_str})"))

    # Make some new files to hold the evaluation results:
    base_results_filename = prediction_file.rstrip(".txt")
    results_path = f"{results_folder}/{base_results_filename}_results.txt"
    precision_recall_path = f"{results_folder}/precision_recall/{base_results_filename}_precision_recall.txt"

    with open(results_path, "w") as results_handle, open(precision_recall_path, "w") as precision_recall_handle:

        results_handle.write(f"AUTHOR/TEAM: {team}")
        results_handle.write(f"MODEL: {model}")
        results_handle.write(f"KEYWORDS: {keywords}")
        results_handle.write(f"Species: {taxonomy}")
        results_handle.write("Ontology\tType\tMode\t | Fmax\tThreshold\tCoverage\n")

        indent = " " * 3

        for ontology in ("bpo", "cco", "mfo"):
            ontology_path = f"{results_folder}/{team}_{model}_{taxonomy}_{ontology.upper()}.txt"
            print(cf.blue("============================\n"))
            print("ONTOLOGY:", cf.green(ontology.upper()))
            for benchmark_type in ("type1", "type2"):

                benchmark_type_str = typeConverter(benchmark_type)

                print(f"{indent}BENCHMARK TYPE:", cf.green(f"{benchmark_type_str} ({benchmark_type})"))
                benchmark, obo_count_dict = read_benchmark(
                    ontology,
                    taxonomy_str,
                    benchmark_type,
                    benchmark_folder,
                    obo_path
                )
                if benchmark is None:
                    sys.stderr.write("No benchmark is available for the input species and type")

                c = PrecREC(benchmark, ontology_path, obo_count_dict[ontology])

                if c.exist:
                    for mode in ("partial", "full"):
                        print(indent*2, "MODE:", cf.green(mode))
                        precision, recall, fmax, threshold = c.Fmax_output(mode)
                        coverage = c.coverage()
                        print(indent*3, "      FMAX:", cf.green(fmax))
                        print(indent*3, " THRESHOLD:", cf.green(threshold))
                        print(indent*3, "  COVERAGE:", cf.green(coverage))
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
