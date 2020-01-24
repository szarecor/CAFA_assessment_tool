from collections import defaultdict
import os
import numpy

valid_species = [
    "all",
    "eukarya",
    "prokarya",
    "HELPY",
    "ECOLI",
    "RAT",
    "DANRE",
    "SULSO",
    "DROME",
    "PSEPK",
    "STRPN",
    "PSEAE",
    "BACSU",
    "MYCGE",
    "HUMAN",
    "METJA",
    "DICDI",
    "YEAST",
    "SCHPO",
    "ARATH",
    "XENLA",
    "MOUSE",
    "PSESM",
    "SALTY",
    "CANAX",
    "SALCH",
]
valid_types = ["type1", "type2", "all"]
valid_subtypes = ["easy", "hard"]
# easy and hard are only in NK benchmarks!!!!
valid_modes = ["full", "partial"]

root_terms = (
    "GO:0008150",  # biological process
    "GO:0005575",  # cellular component
    "GO:0003674",  # molecular function
)


class PrecisionRecallCalculator:
    """ A class for doing precision recall calculations """

    def __init__(self, benchmark, ontology_specific_prediction_path, ontology_term_count):
        """
        Initalize.

        Input:
        benchmark -- instance of the benchmark class
        ontology_specific_prediction_path_path -- !ontology-specific! prediction file separated in GOPred (without headers)
        obocounts -- total number of terms in the ontology
        """
        self.exists = True
        self.ancestors = benchmark.ancestors
        self.true_terms = benchmark.true_terms
        self.ontology_term_count = ontology_term_count
        # flag corresponding if the program has been run, allowing print functions
        self.ran = False
        # Set of all obsolete terms found
        self.obsolete_terms = set()
        # count_above_threshold is the number of proteins with at least one term above threshold
        self.count_above_threshold = defaultdict()
        # count_predictions_in_benchmark is the number of predicted proteins in this file that are in the benchmark file (for coverage)
        self.count_predictions_in_benchmark = 0
        # self.data is the dictionary for all predicted terms
        self.data = defaultdict(list)
        # self.predicted_bench is the dictionary for all predicted terms that are benchmarks
        self.predicted_bench = defaultdict(defaultdict)
        # key:protein
        # value: list of dictionaries
        # key: GO term
        # value: tuple(confidence, Boolean) Boolean is whether in self.true_terms

        # Read in prediction file

        #print("WHAT IS PRED PATH?", ontology_specific_prediction_path)


        #if os.path.getsize(ontology_specific_prediction_path) > 0:
        for prediction_line in open(ontology_specific_prediction_path, "r"):

            target, go_id, confidence = [i.strip() for i in prediction_line.split()]
            self.data[target].append(
                {"term": go_id, "confidence": float(confidence)}
            )

        # Propagated prediction
        for protein in self.data:

            '''
            print("!!!!", protein, self.true_terms[protein])

            print("TRUE TERMS:")
            print(self.true_terms)
            print("END TRUE TERMS")
            print(protein in self.true_terms.keys(), self.true_terms[protein])
            print(self.data[protein])
            print("==========================")
            '''


            if self.true_terms[protein]:
                """
                benchmark.true_terms[protein] not an empty set
                The protein is in the benchmark file
                i.e. gained experimental annota
                """
                self.count_predictions_in_benchmark += 1
                #print("incrementing", self.count_predictions_in_benchmark)

                '''
                import pprint
                pprint.pprint(self.data[protein])
                print("=" * 20)
                '''

                # Renaming var here, I _think_ tc is term_confidence:
                for term_confidence in self.data[protein]:
                # for tc in self.data[protein]:

                    term = term_confidence.get("term")
                    confidence = term_confidence.get("confidence")

                    try:
                        ancestor_terms = self.ancestors[term].difference(root_terms)
                    except KeyError:
                        # Add unknown term to the obsolete set
                        self.obsolete_terms.add(term)
                        continue
                    # For each term
                    if term in self.predicted_bench[protein]:
                        # Term already exists, update confidence
                        self.update_confidence(protein, term_confidence)
                    else:
                        # Add term to self.predicted_bench
                        # Add confidence and compare with self.true_terms
                        # Regardless of comparision, propagate
                        if term not in root_terms:
                            self.predicted_bench[protein][term] = self.is_valid_term(protein, term_confidence)
                            for ancestor_term in ancestor_terms:
                                # Make a new TC with ancestors
                                new_term_confidence = {
                                    "term": ancestor_term,
                                    "confidence": confidence,
                                }
                                if ancestor_term in self.predicted_bench[protein]:
                                    # Term already exists, update confidence
                                    self.update_confidence(protein, new_term_confidence)
                                else:
                                    # Add term to self.predicted_bench
                                    self.predicted_bench[protein][ancestor_term] = self.is_valid_term(protein, new_term_confidence)

            if self.count_predictions_in_benchmark == 0:
                self.exists = False
                '''
                #print(self.data)
                print("TRUE TERMS:")
                print(self.true_terms)
                print("TRUE TERM COUNT:")
                print(self.true_term_count)
                print("OBSOLETE TERMS:")
                print(self.obsolete_terms)
                '''

                print("No protein in this predicted set became a benchmark\n")
        #else:
        #    self.exists = False
        #    print("No prediction made in this ontology.\n")

    @property
    def true_term_count(self):
        return len(self.true_terms)

    @property
    def coverage(self):
        """
        Determine the coverage.

        Finds the value of the coverage-> BETTER EXPLANATION?

        Output:
        The coverage as a float
        """
        return float(self.count_predictions_in_benchmark) / self.true_term_count

    def update_confidence(self, protein: str, term_confidence: dict):
        """
        Update Confidence for given protein and propagate.

        This function compares the confidence value in tc to the confidence in self.predicted
        If tc is larger, than it overwrites the confidence in self.predicted

        Input:
        protein -- chosen protein
        tc -- a dictionary format as: {'confidence':0.57,'term':'GO:006644'}

        Output:
        None
        """
        # Defined for readablity
        term = term_confidence.get("term")
        confidence = term_confidence.get("confidence")

        #print("====", self.predicted_bench[protein][term])  # [term][0] = confidence

        if confidence > self.predicted_bench[protein][term][0]:
            # Update the confidence
            #print("================================")

            # Not sure what's going on here.
            # I _might_ have broken something converting a list to a tuple elsewhere in the code...
            #self.predicted_bench[protein][term][0] = confidence
            self.predicted_bench[protein][term] = (confidence, self.predicted_bench[protein][term][1])
            # Propagate changes if necessary
            for ancestor_term in self.ancestors[term].difference(root_terms):
                if confidence > self.predicted_bench[protein][ancestor_term][0]:
                    # Update the confidence
                    #self.predicted_bench[protein][ancestor_term][0] = confidence
                    self.predicted_bench[protein][ancestor_term] = (confidence, self.predicted_bench[protein][term][1])

    def is_valid_term(self, protein: str, term_confidence: dict) -> iter:
        """
        Check if tc['term'] is a True term.

        This function compares if tc['term'] is in self.true_terms

        Input:
        protein -- chosen protein to for this comparision
        tc -- dictionary with {'confidence':0.57,'term':'GO:006644'}

        Output:
        A list containing:
        'confidence'
        Boolean Value
        """

        '''
        if term_confidence.get("term") in self.true_terms.get(protein):
            return (term_confidence.get("confidence"), True)
        else:
            return (term_confidence.get("confidence"), False)
        '''
        is_valid = term_confidence.get("term") in self.true_terms.get(protein)
        return term_confidence.get("confidence"), is_valid

    '''
   @property
   def obsolete_terms(self):
      """ Get obsolete terms used by the prediction team. """
      return self.obsolete_terms
   '''

    def term_precision_recall(self, threshold, protein):
        """
         Calculate the precision recall of a single protein.

         Input:
         threshold -- the value for PRRC threshold
         protein -- the chosen protein
      """

        # Initalize Variables
        true_positive_count = 0  # True positive
        count = 0  # Count how many terms are above the threshold

        if threshold == 0:
            # At threshold 0, every term (propagated) in the ontology is considered predicted
            # so TP is all the true terms in the benchmark
            # recall = 1
            # precision is TP over all terms in the ontology
            true_positive_count = float(len(self.true_terms[protein]))
            count = self.ontology_term_count
        else:
            for term in self.predicted_bench[protein]:
                # If it is above the threshold, increment the count
                if self.predicted_bench[protein][term][0] >= threshold:
                    count += 1

                    # If it is actually True, increment TP
                    if self.predicted_bench[protein][term][1]:
                        true_positive_count += 1
        # Find PR: TP / (TP + FP)
        try:
            precision = true_positive_count / count
        except ZeroDivisionError:
            precision = None

        # Find RC: TP (TP + FN)
        recall = true_positive_count / len(self.true_terms[protein])
        # Safe becuase protein is in the benchmark file
        #print(precision, recall, count, true_positive_count)
        return precision, recall, count, true_positive_count

    def precision_recall(self, threshold, mode):
        """
      Calculate the overall (average) PRRC of the team.

      Input:
      threshold -- The threshold for PRRC
      'partial':  recall is averaged over all benchmark proteins that are predicted

      Output:
      precision --
      recall --
      """

        # Initialize Variables
        PR = 0.0
        RC = 0.0
        self.count_above_threshold[threshold] = 0

        for protein in self.predicted_bench:
            pr, rc = self.term_precision_recall(threshold, protein)[0:2]

            #print("what is result?", pr, rc)

            if pr is not None:
                PR += pr
                self.count_above_threshold[threshold] += 1
            if rc is not None:
                RC += rc

        if mode == "partial":
            try:
                recall = RC / self.count_predictions_in_benchmark
            except ZeroDivisionError:
                recall = 0
                print("No protein in this predicted set became benchmarks\n")

        elif mode == "full":
            try:
                recall = RC / self.true_term_count  # count_true_terms
            except ZeroDivisionError:
                recall = 0
                print("No protein in this benchmark set\n")

        try:
            precision = PR / self.count_above_threshold[threshold]
        except ZeroDivisionError:
            precision = None
            print("No prediction is made above the %.2f threshold\n" % threshold)

        # PRRC has run
        self.ran = True
        return (precision, recall)

    def Fmax_output(self, mode):
        """
      Compute PRRC for every threshold and Fmax value.

      Input:
      mode -- The evaluation mode used ('full' or 'partial')

      Output:
      A list containing:
      PR              -- A list containing the precison values
      RC              -- A list containing the recall values
      fmax            -- the maximun f value found
      fmax_threshold  -- The threshold that corresponds to fmax
      """

        # Intialize Variables
        fmax = 0.0
        fmax_threshold = 0.0
        PR = []
        RC = []

        # Run over all threshold values from 0 to 1, two signifigant digits
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):

            threshold = numpy.around(threshold, decimals=2)
            # Run PRRC on given threshold
            pr, rc = self.precision_recall(threshold, mode)
            if pr is None:
                # No prediction above this threshold
                break
            else:
                PR.append(pr)
                RC.append(rc)
                # Find the F-value for this particular threshold
                try:
                    f = (2 * pr * rc) / (pr + rc)
                except ZeroDivisionError:
                    f = None

            if f is not None and f >= fmax:  ###########QUESTION##############
                fmax = f
                fmax_threshold = threshold
        # Have found the Fmax at this point
        return [PR, RC, fmax, fmax_threshold]

    '''
   def print_protein_count(self, threshold):
      """
      Print to console the counts if PRRC has been run.

      Input:
      threshold --The threshold value for count_above_threshold
      """

      if self.ran is True:
         print("number of benchmark proteins: %s\n" % self.count_true_terms)
         # Those with not-None recall:
         # (predicted protein that are in benchmark, only one species per prediction file!! )
         print(
            "number of predicted proteins that are in the benchmark file: %s\n"
            % self.count_predictions_in_benchmark
         )
         # Those with not-None precision:
         print(
            "number of proteins with at least one term above threshold: %s\n"
            % self.count_above_threshold[threshold]
         )

      else:
         print("Run precision_recall(%s) first\n" % str(threshold))

   def print_confidence(self, output_path):
      """
      Print confidence and True/False to a file.

      Input:
      output_path -- Where to save the file
      """

      protindex = 0
      out = open(output_path, "w")
      for prot in self.predicted_bench:
         protindex += 1
         for term in self.predicted_bench[prot]:
            out.write(
               "%s\t%s\t%s\n"
               % (
                  str(protindex),
                  self.predicted_bench[prot][term][0],
                  self.predicted_bench[prot][term][1],
               )
            )
      out.close()
   '''
