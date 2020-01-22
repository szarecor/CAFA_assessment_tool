#!/usr/bin/env python4

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

import re
import os
from dataclasses import dataclass
from typing import IO
from collections import defaultdict
from Ontology.IO import OboIO



@dataclass
class Prediction:
    target: str = None
    go_term: str = None
    confidence: float = None


pr_field = re.compile("^PR=[0,1]\.[0-9][0-9];$")
rc_field = re.compile("^RC=[0,1]\.[0-9][0-9]$")
# Fix to add HP ontology 2014-1-9
go_field = re.compile("^(GO|HP):[0-9]{5,7}$")
# Fix to add EFI targets 2014-1-9
target_field = re.compile("^(T|EFI)[0-9]{5,20}$")
confidence_field = re.compile("^[0,1]\.[0-9][0-9]$")

# Legal states: the CAFA prediction records fields, and their order. KEYWORDS and ACCURACY are
# optional
legal_states1 = ["author", "model", "keywords",
                 "accuracy", "go_prediction", "end"]
legal_states2 = ["author", "model", "keywords", "go_prediction", "end"]
legal_states3 = ["author", "model", "go_prediction", "end"]
legal_states4 = ["go_prediction"]

legal_keywords = [
    "sequence alignment", "sequence-profile alignment", "profile-profile alignment", "phylogeny",
    "sequence properties",
    "physicochemical properties", "predicted properties", "protein interactions", "gene expression",
    "mass spectrometry",
    "genetic interactions", "protein structure", "literature", "genomic context", "synteny",
    "structure alignment",
    "comparative model", "predicted protein structure", "de novo prediction", "machine learning",
    "genome environment",
    "operon", "ortholog", "paralog", "homolog", "hidden markov model", "clinical data", "genetic data",
    "natural language processing", "other functional information"
]


class GOPrediction:
    """
    A class for reading and storing CAFA GO predictions
    self.data is a dictionary
       key: protein ID
       value: [{'term':go_term_1, 'threshold': threshold_1},...,{'term':go_term_n, 'threshold': threshold_n}]
    in self.read(pred_path) , pred_path should be a handle
    Edited by Ashley: 05/26/2016
    The read function should read in from a RAW prediction file submitted by a CAFA participanting team (post format check though)
    The split function will perform spliting the prediction by ontology, which used to be a separate function in preprocess.py
    updated by Ashley: 12/29/2016
    """

    def __init__(self):
        # get author (teamID) and model from both header and filename
        self.author = None
        self.model = None
        # get keywords and accuracy from header
        self.keywords = []
        # split by ontology
        # get taxon from filename
        self.taxon = None
        self.data = defaultdict(list)

    def _author_check(self, inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 2:
            correct = False
            errmsg = "AUTHOR: invalid number of fields. Should be 2"
        elif fields[0] != "AUTHOR":
            correct = False
            errmsg = "AUTHOR: First field should be AUTHOR"
        else:
            self.author = fields[1]
        return correct, errmsg

    def _model_check(self, inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 2:
            correct = False
            errmsg = "MODEL: invalid number of fields. Should be 2"
        elif fields[0] != "MODEL":
            correct = False
            errmsg = "MODEL: First field should be MODEL"
        elif len(fields[1]) != 1 or not fields[1].isdigit():
            correct = False
            errmsg = "MODEL: second field should be single digit."
        else:
            self.model = int(fields[1])
        return correct, errmsg

    def _keywords_check(self, inrec):
        correct = True
        errmsg = None
        if inrec[:8] != "KEYWORDS":
            correct = False
            errmsg = "KEYWORDS: first field should be KEYWORDS"
        else:
            keywords = [i.strip().lower() for i in inrec[8:].split(",")]
            for keyword in keywords:
                # stupid full stop
                if keyword[-1] == ".":
                    keyword = keyword[:-1]
                if keyword not in legal_keywords:
                    correct = False
                    errmsg = "KEYWORDS: illegal keyword %s" % keyword
                    break
                else:
                    self.keywords.append(keyword)
        return correct, errmsg

    def _accuracy_check(self, inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 4:
            correct = False
            errmsg = "ACCURACY: error in number of fields. Should be 4"
        elif fields[0] != "ACCURACY":
            correct = False
            errmsg = "ACCURACY: first field should be 'ACCURACY'"
        elif not fields[1].isdigit() or len(fields[1]) != 1:
            correct = False
            errmsg = "ACCURACY: second field should be a single digit"
        elif not pr_field.match(fields[2]):
            correct = False
            errmsg = "ACCURACY: error in PR field"
        elif not rc_field.match(fields[3]):
            correct = False
            errmsg = "ACCURACY: error in RC field"
        return correct, errmsg

    def _go_prediction_check(self, inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 3:
            correct = False
            errmsg = "GO prediction: wrong number of fields. Should be 3"
        elif not target_field.match(fields[0]):
            correct = False
            errmsg = "GO prediction: error in first (Target ID) field"
        elif not go_field.match(fields[1]):
            correct = False
            errmsg = "GO prediction: error in second (GO ID) field"
        elif not confidence_field.match(fields[2]):
            correct = False
            errmsg = "GO prediction: error in third (confidence) field"
        elif float(fields[2]) > 1.0:
            correct = False
            errmsg = "GO prediction: error in third (confidence) field. Cannot be > 1.0"
        else:
            self.data[fields[0]].append(
                {'term': fields[1], 'confidence': float(fields[2])})
        return correct, errmsg

    def _end_check(self, inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 1:
            correct = False
            errmsg = "END: wrong number of fields. Should be 1"
        elif fields[0] != "END":
            correct = False
            errmsg = "END: record should include the word END only"
        return correct, errmsg

    def _handle_error(self, correct, errmsg, inrec):
        if not correct:
            print(inrec)
            raise ValueError(errmsg)

    def parse_prediction_file(self, prediction_handle: IO):
        prediction_handle.seek(0)

        for line in prediction_handle:
            line_data = self.parse_prediction_line(line)
            if line_data is not None:
                self.data[line_data.target].append(line_data)

        return self.data

    def parse_prediction_line(self, line_str: str) -> Prediction:
        '''
        Parses a line of CAFA prediction text and returns a named tuple
        containing target, go term ID and confidence
        Input expected in this form:
        T96060004121 GO:0008152 0.01

        If input does not match expectations, this function returns the None value
        '''

        if type(line_str) is bytes:
            line_str = line_str.decode()

        try:
            target_id, go_id, confidence = line_str.strip().split()
        except ValueError:
            # Wrong number of things in the split string!
            return None

        try:
            confidence = float(confidence)
        except ValueError:
            # confidence isn't a number...
            return None

        return Prediction(target_id, go_id, float(confidence))


    def read(self, pred_path):

        visited_states = []
        n_accuracy = 0
        first_prediction = True
        first_accuracy = True
        first_keywords = True
        n_models = 0
        filename = pred_path.name.split('/')[-1]
        filenamefields = filename.split('.')[0].split('_')
        self.taxon = int(filenamefields[2])
        #inline = open(pred_path)
        for inline in pred_path:
            # gzipped files are in bytes. Need to convert to utf-8
            if type(inline) is bytes:
                inline = inline.decode("utf-8")
            inrec = [i.strip() for i in inline.split()]


            field1 = inrec[0]
            # Check which field type (state) we are in
            if field1 == "AUTHOR":
                state = "author"
            elif field1 == "MODEL":
                state = "model"
            elif field1 == "KEYWORDS":
                state = "keywords"
            elif field1 == "ACCURACY":
                state = "accuracy"
            elif field1 == "END":
                state = "end"
            else:  # default to prediction state
                state = "go_prediction"
            # Check for errors according to state
            if state == "author":
                correct, errmsg = self._author_check(inline)
                self._handle_error(correct, errmsg, inline)
                visited_states.append(state)

            elif state == "model":
                n_models += 1
                n_accuracy = 0
                if n_models > 3:
                    raise ValueError("Too many models. Only up to 3 allowed")
                correct, errmsg = self._model_check(inline)
                if correct and self.model != int(filenamefields[1]):
                    correct = False
                    errmsg = 'MODEL: model number in file different from filename: %' % pred_path.name
                self._handle_error(correct, errmsg, inline)
                if n_models == 1:
                    visited_states.append(state)
            elif state == "keywords":
                if first_keywords:
                    visited_states.append(state)
                    first_keywords = False
                correct, errmsg = self._keywords_check(inline)
                self._handle_error(correct, errmsg, inline)
            elif state == "accuracy":
                if first_accuracy:
                    visited_states.append(state)
                    first_accuracy = False
                n_accuracy += 1
                if n_accuracy > 3:
                    self._handle_error(
                        False, "ACCURACY: too many ACCURACY records: %" % pred_path.name)
                else:
                    correct, errmsg = self._accuracy_check(inline)
            elif state == "go_prediction":
                correct, errmsg = self._go_prediction_check(inline)
                self._handle_error(correct, errmsg, inline)
                if first_prediction:
                    visited_states.append(state)
                    first_prediction = False
            elif state == "end":
                correct, errmsg = self._end_check(inline)
                self._handle_error(correct, errmsg, inline)
                visited_states.append(state)
            # End file forloop
        if (visited_states != legal_states1 and
            visited_states != legal_states2 and
            visited_states != legal_states3 and
                visited_states != legal_states4):
            print(visited_states)
            print("file not formatted according to CAFA specs")
            print("Check whether all these record types are in your file")
            print(
                "Check whether all these record types are in your file in the correct order")
            print("AUTHOR, MODEL, KEYWORDS, ACCURACY (optional), predictions, END")
            raise ValueError



    def split_predictions_by_namespace(self, obo_path, prediction_handle):

        molecular_function = []
        biological_process = []
        cellular_component = []

        obo_handle = open(obo_path, "r")
        ontology_reader = OboIO.OboReader(obo_handle).read()

        if self.data is None or len(self.data.keys()) == 0:
            self.parse_prediction_file(prediction_handle)

        for k, predictions in self.data.items():

            for prediction in predictions:

                namespace = ontology_reader.get_namespace(prediction.go_term)

                if namespace == 'molecular_function':
                    molecular_function.append(prediction)
                elif namespace == 'biological_process':
                    biological_process.append(prediction) #= biological_process + v
                elif namespace == 'cellular_component':
                    cellular_component.append(prediction) #= cellular_component + v

        return {
            'molecular function': molecular_function,
            'biological process': biological_process,
            'cellular component': cellular_component,
        }


    def read_and_split_and_write(self, obo_path, pred_path):
        # This function has self.read() included
        # pred_path should be a handle
        # split both writes to the predictions to three separate files
        # also saves the data in memory?
        self.read(pred_path)
        go_graph = OboIO.OboReader(open(obo_path)).read()
        mfo_out = open("%s_MFO.txt" % os.path.splitext(pred_path.name)[0], "w")
        bpo_out = open("%s_BPO.txt" % os.path.splitext(pred_path.name)[0], "w")
        cco_out = open("%s_CCO.txt" % os.path.splitext(pred_path.name)[0], "w")
        for protein in self.data.items():
            for u in protein[1]:
                try:
                    namespace = go_graph.get_namespace(u['term'])
                except KeyError:
                    #print("Term %s not included in obo file.\n" % str(e))
                    namespace = None
                    # This excludes a term that's not found in the obo file
                    # to be entered into any of the ontology-specific files
                if namespace == 'molecular_function':
                    mfo_out.write("%s\t%s\t%.2f\n" %
                                  (protein[0], u['term'], u['confidence']))
                elif namespace == 'biological_process':
                    bpo_out.write("%s\t%s\t%.2f\n" %
                                  (protein[0], u['term'], u['confidence']))
                elif namespace == 'cellular_component':
                    cco_out.write("%s\t%s\t%.2f\n" %
                                  (protein[0], u['term'], u['confidence']))
                elif namespace == None:
                    #print("Some terms are not included in obo file\n")
                    continue
                else:
                    raise ValueError(
                        "Term %s not found in any ontology" % u['term'])
        mfo_out.close()
        bpo_out.close()
        cco_out.close()
        #self.data = None
        # the ontology-specific files will be re-read from PrecREC
        # 12/31/2016


if __name__ == "__main__":
    import sys
    print(sys.path)
    sys.path.append('../')
    sys.path.append('.')

    test_data_path = "./test/test_data/"
    obo_path = "{}go_cafa3.obo".format(test_data_path)
    prediction_path = "ZZZ_1_9606.txt".format(test_data_path)

    test = GOPrediction()

    test.split_predictions_by_namespace(obo_path, prediction_path)